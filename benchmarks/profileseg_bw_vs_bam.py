"""bigWig vs BAM input on the same simulated fragments (ProfileSeg v3).

For each simulated experiment, fragment midpoints are drawn from the simulator's per-bin counts and
given lengths ~ N(150, 25) bp. The same fragments are written twice: as a paired-end BAM (proper
pairs, 50 bp reads) and as a CPM-scaled fragment-coverage bigWig (as most pipelines would produce).
Both inputs go through `extract` → `analyse` → `per_region_tests` with the BAM path in midpoint mode.

Reported per input: null and amplitude-only false-call shares, shift / broadening power, FDP, agreement
of the per-region calls between inputs, and the depth estimate of the bigWig path vs the exact BAM depth.

    python benchmarks/profileseg_bw_vs_bam.py --out DIR [--reps 4]
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config, extract
from genomecanvas.profileseg.regional import per_region_tests
from genomecanvas.profileseg.simulate import SimSettings, simulate

SCENARIOS = ("null", "amplitude", "site_amplitude", "shift", "broaden")
SPACING = 10_000
READ = 50


def write_inputs(prof, lib, out: Path, rng) -> list[dict]:
    import pyBigWig
    import pysam
    edges = prof.edges
    width = np.diff(edges)
    anchors = np.array([r["start"] for r in prof.regions])
    chrom_len = int(anchors.max() + SPACING)
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chrSim", "LN": chrom_len}]}
    samples = []
    for i, smp in enumerate(prof.samples):
        counts = np.rint(np.nan_to_num(prof.values[i]) * lib[i]).astype(int)
        mids, lens = [], []
        for r in range(counts.shape[0]):
            nz = np.flatnonzero(counts[r])
            k = counts[r, nz]
            m = anchors[r] + np.repeat(edges[:-1][nz], k) + rng.uniform(0, 1, k.sum()) * np.repeat(width[nz], k)
            mids.append(m)
        mids = np.concatenate(mids).astype(int)
        lens = np.clip(rng.normal(150, 25, len(mids)).astype(int), 60, 400)
        starts = mids - lens // 2
        ok = starts > 0
        mids, lens, starts = mids[ok], lens[ok], starts[ok]
        # BAM: proper pairs
        unsorted = out / f"{smp['sample_id']}.unsorted.bam"
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
            for j, (s, L) in enumerate(zip(starts.tolist(), lens.tolist(), strict=True)):
                for first in (True, False):
                    a = pysam.AlignedSegment()
                    a.query_name = f"f{j}"
                    a.query_sequence = "A" * READ
                    a.flag = (99 if first else 147)
                    a.reference_id = 0
                    a.reference_start = s if first else s + L - READ
                    a.mapping_quality = 60
                    a.cigarstring = f"{READ}M"
                    a.next_reference_id = 0
                    a.next_reference_start = s + L - READ if first else s
                    a.template_length = L if first else -L
                    a.query_qualities = pysam.qualitystring_to_array("I" * READ)
                    bam.write(a)
        bam_path = out / f"{smp['sample_id']}.bam"
        pysam.sort("-o", str(bam_path), str(unsorted))
        pysam.index(str(bam_path))
        unsorted.unlink()
        # bigWig: CPM-scaled fragment coverage
        cov = np.zeros(chrom_len + 1)
        np.add.at(cov, starts, 1.0)
        np.add.at(cov, starts + lens, -1.0)
        cov = np.cumsum(cov)[:chrom_len] * (1e6 / len(starts))
        change = np.flatnonzero(np.diff(cov)) + 1
        bounds = np.r_[0, change, chrom_len]
        vals = cov[bounds[:-1]]
        keep = vals > 0
        bw_path = out / f"{smp['sample_id']}.bw"
        bw = pyBigWig.open(str(bw_path), "w")
        bw.addHeader([("chrSim", chrom_len)])
        bw.addEntries(["chrSim"] * int(keep.sum()), bounds[:-1][keep].tolist(), ends=bounds[1:][keep].tolist(),
                      values=vals[keep].tolist())
        bw.close()
        samples.append(dict(sample_id=smp["sample_id"], condition=smp["condition"], bam=str(bam_path), bw=str(bw_path)))
    return samples


def run(job):
    scenario, n, rep, outdir = job
    rng = np.random.default_rng(8_000_000 + 1000 * SCENARIOS.index(scenario) + 100 * n + rep)
    st = SimSettings(n_regions=200, n_ref=n, n_alt=n)
    prof, truth = simulate(scenario, st, seed=int(rng.integers(1 << 30)))
    for k, r in enumerate(prof.regions):
        r.update(chrom="chrSim", start=SPACING * (k + 1), end=SPACING * (k + 1) + 1, strand="+")
    work = Path(tempfile.mkdtemp(prefix=f"{scenario}_{n}_{rep}_", dir=outdir))
    samples = write_inputs(prof, np.asarray(truth["library"], float), work, rng)
    changed = (np.asarray(truth["regions"]["shift"]) != 0) | (scenario == "broaden")
    cfg = Config(reference="WT", alternative="KO", anchor="start", orientation="oriented", flank=st.flank,
                 binsize=st.binsize, background_width=500, core_flank=1000)
    res, calls = {}, {}
    for kind in ("bam", "bw"):
        smp = [dict(sample_id=s["sample_id"], condition=s["condition"], signal_file=s[kind]) for s in samples]
        p = extract(cfg, smp, prof.regions)
        reg = per_region_tests(analyse(p, cfg, bootstrap=0))
        out = {}
        for m in ("shape_kernel", "shift"):
            t = reg.tests[m]
            c = np.isin(t["call"], ("Higher", "Lower", "Changed"))
            tr = changed[t["region_index"]]
            nc = int(c.sum())
            out[m] = dict(called=float(c.mean()), fdp=float((c & ~tr).sum() / nc) if nc else 0.0,
                          power=float((c & tr).sum() / tr.sum()) if tr.any() else float("nan"))
            calls[(kind, m)] = dict(zip(t["region_index"].tolist(), c.tolist(), strict=True))
        out["depth_median"] = float(np.median(list(reg.depth.get("median_fragments", {}).values()) or [np.nan]))
        res[kind] = out
    for m in ("shape_kernel", "shift"):
        a, b = calls[("bam", m)], calls[("bw", m)]
        common = sorted(set(a) & set(b))
        both = sum(a[k] and b[k] for k in common)
        either = sum(a[k] or b[k] for k in common)
        res[f"agreement_{m}"] = both / either if either else float("nan")
    for f in work.iterdir():
        f.unlink()
    work.rmdir()
    return scenario, n, res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=5)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    jobs = [(sc, n, r, str(out)) for sc in SCENARIOS for n in (2, 3) for r in range(args.reps)]
    agg = {}
    with ProcessPoolExecutor(args.jobs) as ex:
        for sc, n, res in ex.map(run, jobs):
            agg.setdefault((sc, n), []).append(res)
    lines = ["# bigWig vs BAM on the same simulated fragments", "",
             "| scenario | design | input | kernel: called / FDP / power | shift: called / FDP / power | median fragments per region (depth diagnostic) |",
             "|---|---|---|---|---|---|"]
    table = {}
    for (sc, n), rows in sorted(agg.items()):
        for kind in ("bam", "bw"):
            m = lambda t, k: float(np.nanmean([r[kind][t][k] for r in rows]))  # noqa: E731
            dep = float(np.nanmean([r[kind]["depth_median"] for r in rows]))
            lines.append(f"| {sc} | {n}v{n} | {kind} | {m('shape_kernel', 'called'):.4f} / {m('shape_kernel', 'fdp'):.3f} / "
                         f"{m('shape_kernel', 'power'):.2f} | {m('shift', 'called'):.4f} / {m('shift', 'fdp'):.3f} / {m('shift', 'power'):.2f} | {dep:.0f} |")
            table[f"{sc}:{n}:{kind}"] = rows
        ag = {m: float(np.nanmean([r[f"agreement_{m}"] for r in rows])) for m in ("shape_kernel", "shift")}
        lines.append(f"| {sc} | {n}v{n} | agreement (Jaccard of calls) | {ag['shape_kernel']:.2f} | {ag['shift']:.2f} | |")
    (out / "bw_vs_bam.json").write_text(json.dumps({k: v for k, v in table.items()}, indent=1, default=float))
    (out / "BW_VS_BAM.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
