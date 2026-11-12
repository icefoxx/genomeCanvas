"""Data survey for agents: what is in a folder, what is wrong with it, and which figures fit.

Used by ``gcplot.py survey`` and the Studio chat agent. Reads only headers, sample sheets and BED
rows (never BigWig values), so it takes seconds even on multi-GB collections. Recommendations are
ready-to-run requests with the reason for each and the questions an agent should put to the user.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

# chr1 (or the first-listed autosome) length -> assembly
GENOMES = {
    248956422: "hg38",
    249250621: "hg19",
    195471971: "mm10",
    195154279: "mm39",
    260522016: "rn7",
    282763074: "rn6",
    23513712: "dm6",
    23011544: "dm3",
    15072434: "ce11",
    15072423: "ce10",
    59578282: "danRer11",
    230218: "sacCer3",
    30427671: "TAIR10",
}
FIRST_CHROM = ("chr1", "1", "chr2L", "2L", "chrI", "I", "Chr1")
MARK = re.compile(
    r"(H[234][AB]?K\d+(?:me[123]|ac|ub\d?)|H2A\.?Z|H3\.3|CTCF|RAD21|SMC[13]|POLR2A|Pol2|Pol_?II|"
    r"RNAPII|ATAC|DNase|MNase|CUT&?(?:RUN|Tag)|input|IgG)",
    re.IGNORECASE,
)
CONTROL = re.compile(r"(^|[_\-.])(input|igg|control|ctrl|wce)([_\-.]|$)", re.IGNORECASE)
REPLICATE = re.compile(r"(?:^|[_\-.])(?:rep|r|replicate)[_\-]?(\d+)(?=[_\-.]|$)", re.IGNORECASE)
GSM = re.compile(r"(GSM\d+)")
MANIFESTS = ("samples.tsv", "genomecanvas.tsv", "manifest.tsv")
MAX_REGIONS = 20000
PARALLEL_SECONDS_PER_GB = 7.5  # measured: 6 mm10 BigWigs, 3.7 GB, 10 kb bins, 4 processes


def survey(paths, suffixes=None, write_dir=None):
    if suffixes is None:
        from .engine import SUFFIXES as suffixes
    files, skipped, folders = [], [], []
    for raw in paths:
        path = Path(raw).expanduser().resolve()
        if path.is_dir():
            folders.append(path)
            for item in sorted(path.rglob("*")):
                rel = item.relative_to(path).parts
                if item.is_file() and not any(p.startswith(".") for p in rel):
                    _classify(item, suffixes, files, skipped)
        elif path.is_file():
            _classify(path, suffixes, files, skipped)
        else:
            skipped.append(dict(path=str(path), reason="not found"))
    project = _project(folders)
    sheet = _sheet(folders)
    for f in files:
        _annotate(f, sheet)
    order = list((sheet or {}).get("rows", {}))
    files.sort(key=lambda f: order.index(f["path"].name) if f["path"].name in order else len(order))
    signals = [f for f in files if f["role"] == "signal"]
    if not sheet:  # no sample sheet: ChIP first, controls last
        signals.sort(key=lambda f: f.get("control", False))
    regions = [f for f in files if f["role"] == "regions"]
    genes = [f for f in files if f["role"] == "genes"]
    checks = _checks(files, signals, regions, project)
    recs = _recommend(signals, regions, genes, project, checks)
    if write_dir:
        _write(recs, Path(write_dir).expanduser().resolve())
    gb = sum(f["size"] for f in signals) / 1024**3
    return dict(
        status="ok",
        project=project,
        sample_sheet=sheet["path"] if sheet else None,
        genome=checks.pop("genome"),
        summary=dict(
            signals=len(signals),
            regions=len(regions),
            genes=len(genes),
            tables=sum(f["role"] == "table" for f in files),
            signal_gb=round(gb, 2),
            marks=sorted({f.get("mark") or f["label"] for f in signals}),
            controls=[f["label"] for f in signals if f.get("control")],
            groups=sorted({f["group"] for f in signals if f.get("group")}),
        ),
        files=[_public(f) for f in files],
        skipped=skipped,
        warnings=checks["warnings"],
        recommendations=recs,
    )


# ----------------------------------------------------------------------------- files


def _classify(path, suffixes, files, skipped):
    if path.name in MANIFESTS or path.name == "project.yaml":
        return
    lower = [s.lower() for s in path.suffixes]
    if lower and lower[-1] == ".gz":
        skipped.append(dict(path=str(path), reason="compressed; decompress it to use it"))
        return
    role = suffixes.get(lower[-1]) if lower else None
    if not role:
        if path.suffix.lower() not in (".txt", ".md", ".log", ".json", ".yaml", ".yml"):
            skipped.append(dict(path=str(path), reason="unsupported file type"))
        return
    files.append(dict(path=path, role=role, size=path.stat().st_size))


def _annotate(f, sheet):
    path = f["path"]
    stem = path.name.split(".")[0]
    row = (sheet or {}).get("rows", {}).get(path.name) or (sheet or {}).get("rows", {}).get(
        str(path)
    )
    f["label"] = (row or {}).get("label") or stem
    f["group"] = (row or {}).get("group") or (row or {}).get("condition") or ""
    if row and row.get("replicate"):
        f["replicate"] = row["replicate"]
    mark = MARK.search(f"{f['label']} {stem}")
    if mark:
        f["mark"] = mark.group(1)
    f["control"] = bool(CONTROL.search(stem) or CONTROL.search(f["label"]))
    rep = REPLICATE.search(stem)
    if rep and "replicate" not in f:
        f["replicate"] = rep.group(1)
    gsm = GSM.search(stem)
    if gsm:
        f["accession"] = gsm.group(1)
    if f["role"] == "signal":
        _signal_header(f)
    elif f["role"] == "regions":
        _regions(f)


def _signal_header(f):
    if f["path"].suffix.lower() not in (".bw", ".bigwig"):
        f["kind"] = "bedGraph"
        return
    try:
        import pyBigWig

        handle = pyBigWig.open(str(f["path"]))
        f["chrom_sizes"] = dict(handle.chroms())
        handle.close()
        f["kind"] = "BigWig"
    except Exception as exc:  # noqa: BLE001 - reported to the agent as a file problem
        f["error"] = f"cannot open BigWig: {exc}"


def _regions(f):
    rows, chroms, best = 0, set(), None
    if f["size"] > 500 * 1024**2:
        f["note"] = "large region file; rows not counted"
        return
    with f["path"].open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            rows += 1
            chroms.add(cols[0])
            # narrowPeak signalValue (col 7), else BED score (col 5): the strongest peak
            score = cols[6] if len(cols) > 6 else cols[4] if len(cols) > 4 else "0"
            try:
                value = float(score)
            except ValueError:
                value = 0.0
            if best is None or value > best[0]:
                best = (value, cols[0], int(cols[1]), int(cols[2]))
    f["rows"], f["chroms"] = rows, sorted(chroms)
    if best:
        f["top"] = dict(chrom=best[1], start=best[2], end=best[3], score=best[0])


def _sheet(folders):
    for folder in folders:
        for name in MANIFESTS:
            path = folder / name
            if not path.is_file():
                continue
            with path.open(encoding="utf-8") as handle:
                lines = [line for line in handle if line.strip() and not line.startswith("#")]
            rows = {}
            for row in csv.DictReader(lines, delimiter="\t"):
                row = {str(k).strip().lower(): (v or "").strip() for k, v in row.items() if k}
                ref = row.get("file") or row.get("path")
                if ref:
                    rows[Path(ref).name] = row
            return dict(path=str(path), rows=rows)
    return None


def _project(folders):
    for folder in folders:
        path = folder / "project.yaml"
        if path.is_file():
            try:
                import yaml

                data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            except Exception as exc:  # noqa: BLE001
                return dict(path=str(path), error=str(exc))
            return dict(
                path=str(path),
                name=data.get("name"),
                genome=data.get("genome"),
                defaults=data.get("defaults") or {},
            )
    return None


def _public(f):
    out = {k: v for k, v in f.items() if k not in ("chrom_sizes", "chroms", "path")}
    out["path"] = str(f["path"])
    if f.get("chrom_sizes"):
        out["n_chroms"] = len(f["chrom_sizes"])
    if f.get("chroms"):
        out["n_chroms"] = len(f["chroms"])
    return out


# ----------------------------------------------------------------------------- checks


def _guess_genome(sizes):
    for name in FIRST_CHROM:
        if name in sizes:
            return GENOMES.get(sizes[name])
    return None


def _checks(files, signals, regions, project):
    warnings = []
    guesses = {}
    for f in signals:
        if f.get("error"):
            warnings.append(f"{f['path'].name}: {f['error']}")
        if f.get("chrom_sizes"):
            guesses.setdefault(_guess_genome(f["chrom_sizes"]) or "unknown", []).append(f["label"])
    known = {g: v for g, v in guesses.items() if g != "unknown"}
    declared = (project or {}).get("genome")
    genome = next(iter(known)) if len(known) == 1 else declared
    if len(known) > 1:
        detail = "; ".join(f"{g}: {', '.join(v)}" for g, v in known.items())
        warnings.append(f"Signals come from different assemblies ({detail}). Do not mix them.")
    if declared and known and declared not in known:
        warnings.append(
            f"project.yaml says {declared} but the BigWig headers look like "
            f"{', '.join(known)}; check before plotting."
        )
    sig_chr = {c.startswith("chr") for f in signals for c in (f.get("chrom_sizes") or {})}
    for f in regions:
        reg_chr = {c.startswith("chr") for c in f.get("chroms", [])}
        if sig_chr and reg_chr and sig_chr != reg_chr:
            warnings.append(
                f"{f['path'].name} names chromosomes differently from the signals "
                "('chr1' vs '1'); GenomeCanvas matches them, but check the assembly."
            )
        if f.get("rows", 0) > MAX_REGIONS:
            warnings.append(
                f"{f['path'].name}: {f['rows']:,} regions; heatmaps keep the {MAX_REGIONS:,} "
                "strongest (narrowPeak signalValue, else BED score)."
            )
    marks = {}
    for f in signals:
        marks.setdefault((f.get("mark") or f["label"]).lower(), []).append(f)
    if signals and all(len(v) == 1 for v in marks.values()):
        warnings.append(
            "One file per mark/factor and no replicates: comparisons between samples are "
            "descriptive only (no statistics)."
        )
    if signals and not any(f.get("control") for f in signals):
        warnings.append("No input/IgG control found; enrichment over background cannot be shown.")
    return dict(warnings=warnings, genome=genome, marks=marks)


# ----------------------------------------------------------------------------- recommend


def _sig(f):
    return dict(path=str(f["path"]), role="signal", label=f["label"], group=f.get("group") or "")


def _reg(f):
    return dict(path=str(f["path"]), role="regions", label=f["label"])


def _locus(regions, project):
    defaults = (project or {}).get("defaults") or {}
    if all(k in defaults for k in ("chrom", "start", "end")):
        return dict(
            chrom=defaults["chrom"],
            start=int(defaults["start"]),
            end=int(defaults["end"]),
            source="project.yaml default locus",
        )
    tops = [f for f in regions if f.get("top")]
    if not tops:
        return None
    top = max(tops, key=lambda f: f["top"]["score"])["top"]
    mid = (top["start"] + top["end"]) // 2
    return dict(
        chrom=top["chrom"],
        start=max(0, mid - 20000),
        end=mid + 20000,
        source=f"strongest peak in {max(tops, key=lambda f: f['top']['score'])['path'].name}",
    )


def _recommend(signals, regions, genes, project, checks):
    recs = []
    chip = [f for f in signals if not f.get("control")]
    controls = [f for f in signals if f.get("control")]
    gb = sum(f["size"] for f in signals) / 1024**3
    title = (project or {}).get("name") or ""
    # The ChIP with its own peak set (e.g. SOX2 + SOX2 peaks) is the natural lead sample.
    peak_names = " ".join(f["label"].lower() for f in regions)
    lead = next((f for f in chip if f["label"].lower() in peak_names), chip[0] if chip else None)

    def add(slug, why, files, settings, ask=(), cost=None, priority=2):
        recs.append(
            dict(
                plot=slug,
                priority=priority,
                why=why,
                ask=list(ask),
                cost=cost,
                request=dict(plot=slug, title=title, files=files, settings=settings),
            )
        )

    if signals and regions:
        anchor = max(regions, key=lambda f: f.get("rows", 0))
        add(
            "enrichment/profile-heatmap",
            f"Signal of every sample around {anchor['label']} ({anchor.get('rows', 0):,} "
            "regions): the standard first "
            "look at where each mark or factor sits relative to the peaks.",
            [_sig(f) for f in signals] + [_reg(anchor)],
            dict(alignment="center", flank=1000, binsize=50),
            ask=[
                f"Anchor on {anchor['label']} (largest set), or on another region set: "
                + ", ".join(f["label"] for f in regions if f is not anchor)
                if len(regions) > 1
                else "Centre on the peaks (default) or use TSS/TES of a gene set?",
                "Leave input/IgG in the figure as a background reference, or drop it?",
            ],
            cost="seconds (matrix cached afterwards)",
            priority=1,
        )
        if len(chip) >= 3:
            add(
                "enrichment/profile-heatmap",
                "Same figure with rows k-means clustered, to find region classes with "
                "different mark combinations (e.g. promoter-like vs enhancer-like).",
                [_sig(f) for f in signals] + [_reg(anchor)],
                dict(
                    alignment="center",
                    flank=1000,
                    binsize=50,
                    cluster="kmeans",
                    clusters=4,
                    # Controls carry no region-specific signal; leave them out of the clustering.
                    cluster_samples=",".join(f["label"] for f in chip),
                    group_label_pos="text-left",
                    group_legend="top",
                ),
                ask=[
                    "How many clusters (k)? 3-6 is usually interpretable.",
                    "Cluster on the ChIP samples only (default here; input/IgG left out) or a "
                    "subset (cluster_samples)?",
                    "After you see the figure: names for the clusters (group_names)?",
                ],
                cost="seconds",
                priority=1,
            )
    if len(signals) >= 2:
        est = gb * PARALLEL_SECONDS_PER_GB
        add(
            "relationships/correlation",
            "Genome-wide similarity between samples (10 kb bins, Spearman): checks that "
            "replicates agree and related marks group together.",
            [_sig(f) for f in signals],
            dict(rel_mode="bins", rel_bin=10000, method="spearman"),
            ask=["Genome-wide bins (default) or only over the peak regions?"],
            cost=f"about {max(5, round(est))} s the first time ({gb:.1f} GB of BigWigs), "
            "then cached",
            priority=1 if len(signals) >= 4 else 2,
        )
        if len(signals) >= 3:
            add(
                "relationships/pca",
                "Same table as the correlation, as a PCA scatter (plus interactive HTML).",
                [_sig(f) for f in signals],
                dict(rel_mode="bins", rel_bin=10000),
                cost="seconds once the correlation table is cached",
            )
    locus = _locus(regions, project)
    if signals and locus:
        add(
            "tracks/stacked",
            f"Browser view of all tracks at one locus ({locus['source']}).",
            [_sig(f) for f in signals] + [_reg(f) for f in regions] + [
                dict(path=str(f["path"]), role="genes", label=f["label"]) for f in genes
            ],
            dict(chrom=locus["chrom"], start=locus["start"], end=locus["end"]),
            ask=["Which gene or locus should the figure show? (suggested: "
                 f"{locus['chrom']}:{locus['start']:,}-{locus['end']:,})"]
            + ([] if genes else ["No GTF found: add one for gene models?"]),
            cost="seconds",
            priority=1 if not regions else 2,
        )
    if chip and controls and locus:
        add(
            "transform/log2-ratio",
            f"{lead['label']} over {controls[0]['label']} at the locus: enrichment over "
            "background.",
            [_sig(lead), _sig(controls[0])],
            dict(chrom=locus["chrom"], start=locus["start"], end=locus["end"]),
            ask=["Which ChIP track over which control?"],
            cost="seconds",
            priority=3,
        )
    for mark, files in checks["marks"].items():
        groups = {}
        for f in files:
            groups.setdefault(f.get("group") or "", []).append(f)
        good = {g: v for g, v in groups.items() if g and len(v) >= 2}
        if len(good) >= 2 and regions:
            (g1, a), (g2, b) = list(good.items())[:2]
            add(
                "differential/volcano",
                f"{mark}: {g1} vs {g2}, {len(a)} vs {len(b)} replicates over the peaks "
                "(exploratory; use DESeq2/DiffBind for formal calls).",
                [_sig(f) for f in a + b] + [_reg(regions[0])],
                dict(numerator=g1, denominator=g2),
                ask=[f"Is {g1} vs {g2} the contrast you want, and which is the reference?"],
                cost="seconds",
                priority=2,
            )
    if lead and len(regions) >= 2:
        add(
            "differential/peak-center",
            f"Compare {lead['label']} at the centres of {regions[0]['label']} vs "
            f"{regions[1]['label']}.",
            [_sig(lead)] + [_reg(f) for f in regions[:2]],
            {},
            ask=["Which signal, and which two region sets?"],
            cost="seconds",
            priority=3,
        )
    if len(signals) >= 12:
        recs.append(
            dict(
                plot="relationships/sample-qc",
                priority=2,
                why=f"{len(signals)} samples: Sample QC flags mislabelled, swapped and "
                "low-quality samples. Needs a sample sheet with mark and group columns and a "
                "bins x samples matrix (gcanvas-sampleqc matrix).",
                ask=["Do you have a sample sheet with mark / group / batch columns?"],
                cost="matrix build: minutes; QC: seconds",
                request=None,
            )
        )
    recs.sort(key=lambda r: r["priority"])
    return recs


def _write(recs, folder):
    import json

    folder.mkdir(parents=True, exist_ok=True)
    for index, rec in enumerate(recs, 1):
        if rec["request"] is None:
            continue
        name = f"{index:02d}-{rec['plot'].replace('/', '-')}.json"
        if rec["request"]["settings"].get("cluster"):
            name = name.replace(".json", "-clustered.json")
        path = folder / name
        path.write_text(json.dumps(rec["request"], indent=2) + "\n", encoding="utf-8")
        rec["request_file"] = str(path)
