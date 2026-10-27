"""Tables for the ProfileSeg review documents, generated from calibration-boundary JSON outputs.

Blocks between ``<!-- gen:NAME -->`` and ``<!-- /gen:NAME -->`` in the target Markdown files are
replaced, so every number in those blocks comes from a recorded run.

    python benchmarks/profileseg_review_tables.py RUNS_DIR docs/PROFILESEG_V3_REVIEW.md [more.md ...]

RUNS_DIR holds one sub-directory per run (dev_v3, dev_v24, dev_fdr_spread, dev_fdr_depth, confirm_v3),
each with calibration_boundary.json.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

TESTS = ("shape_kernel", "shift", "shape_any")
SELECTED = [("width", "60", "2v2"), ("width", "400", "3v3"), ("depth", "5.0", "3v3"), ("noise", "high", "3v3"),
            ("proportion", "0.0", "2v2"), ("proportion", "0.0", "3v3"), ("proportion", "0.5", "2v2"),
            ("correlation", "overlap", "2v2"), ("batch", "balanced 15 bp", "3v3"), ("semireal", "real noise", "2v2"),
            ("mixed_large", "3 classes × 1000", "2v2"), ("mixed_large", "3 classes × 1000", "3v3")]
CLASS_NAMES = {"class0": "narrow (sd 60 bp)", "class1": "broad (sd 400 bp)", "class2": "low depth"}


def load(d: Path) -> dict | None:
    f = d / "calibration_boundary.json"
    return json.loads(f.read_text()) if f.exists() else None


def cell(doc, f, lv, des, label="all"):
    for key, c in doc["cells"].items():
        kf, klv, kd, kl = key.split("|")
        if kf == f and kd == des and kl == label and (klv == lv or _num(klv) == _num(lv)):
            return c
    return None


def _num(x):
    try:
        return float(x)
    except ValueError:
        return x


def fmt(d):
    if "fwer" in d:
        return f"FWER {d['fwer']:.3f} [{d['low']:.3f}, {d['high']:.3f}]"
    return f"{d['fdp']:.3f} [{max(d['low'], 0):.3f}, {d['high']:.3f}] / {d['power']:.2f}"


def block_provenance(runs):
    rows = ["| run | purpose | seeds | null | FDR | commit | uncommitted src files | regional.py sha256 | src diff sha256 | started → finished |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    purpose = dict(dev_v3="development, v3 default", dev_v24="development, v2.4 null", dev_fdr_spread="development, spread FDR strata",
                   dev_fdr_depth="development, depth FDR strata", confirm_v3="confirmatory, v3 default (frozen)")
    for name, doc in runs.items():
        pv, o = doc["provenance"], doc["provenance"]["options"]
        rows.append(f"| {name} | {purpose.get(name, '')} | {o['seed_base']}+ | {o['null']} | {o['fdr_strata']} | `{pv['commit'][:10]}` | "
                    f"{len(pv['uncommitted_src'])} | `{pv['regional_py_sha256']}` | `{pv.get('src_diff_sha256', '–')}` | {doc['started']} → {doc['finished']} |")
    return "\n".join(rows)


def block_wholelist(runs):
    out = ["| run | cells × tests judged | ok | WATCH | EXCEEDS | largest mean FDP | largest CI upper bound |", "|---|---|---|---|---|---|---|"]
    for name in ("dev_v3", "dev_v24", "confirm_v3"):
        doc = runs.get(name)
        if not doc:
            continue
        v, fd, hi = [], [], []
        for key, c in doc["cells"].items():
            if not key.endswith("|all"):
                continue
            for t in TESTS:
                v.append(c[t]["verdict"])
                if "fdp" in c[t]:
                    fd.append(c[t]["fdp"]); hi.append(c[t]["high"])
        out.append(f"| {name} | {len(v)} | {v.count('ok')} | {v.count('WATCH')} | {v.count('EXCEEDS')} | {max(fd):.3f} | {max(hi):.3f} |")
    out += ["", "Selected cells (FDP [95 % CI] / power; whole list):", "",
            "| factor | level | design | run | " + " | ".join(TESTS) + " |", "|---|---|---|---|---|---|---|"]
    for f, lv, des in SELECTED:
        for name in ("dev_v3", "confirm_v3"):
            doc = runs.get(name)
            c = doc and cell(doc, f, lv, des)
            if c:
                out.append(f"| {f} | {lv} | {des} | {name} | " + " | ".join(fmt(c[t]) for t in TESTS) + " |")
    return "\n".join(out)


def block_perclass(runs):
    out = ["Share of **unchanged** regions with p < 0.01 (nominal 0.01), mixed designs:", "",
           "| design | class | test | v2.4 null (dev) | v3 null (dev) | v3 null (confirmatory) |", "|---|---|---|---|---|---|"]
    for f in ("mixed", "mixed_large"):
        for des in ("2v2", "3v3"):
            for label, cname in CLASS_NAMES.items():
                for t in TESTS:
                    vals = []
                    for name in ("dev_v24", "dev_v3", "confirm_v3"):
                        doc = runs.get(name)
                        c = doc and next((cc for k, cc in doc["cells"].items() if k.startswith(f"{f}|") and k.endswith(f"|{des}|{label}")), None)
                        vals.append(f"{c[t]['p01_unchanged']:.4f}" if c else "–")
                    out.append(f"| {f} {des} | {cname} | {t} | " + " | ".join(vals) + " |")
    return "\n".join(out)


def block_fdrstrata(runs):
    out = ["Class FDP (mean over runs) and whole-list FDP [95 % CI], `mixed_large`, by FDR strata:", "",
           "| design | FDR strata | test | narrow | broad | low depth | whole list |", "|---|---|---|---|---|---|---|"]
    for des in ("2v2", "3v3"):
        for name, lab in (("dev_v3", "none (global BH)"), ("dev_fdr_spread", "spread"), ("dev_fdr_depth", "depth")):
            doc = runs.get(name)
            if not doc:
                continue
            for t in TESTS:
                cls = []
                for label in CLASS_NAMES:
                    c = next((cc for k, cc in doc["cells"].items() if k.startswith("mixed_large|") and k.endswith(f"|{des}|{label}")), None)
                    cls.append(f"{c[t]['fdp']:.3f}" if c else "–")
                a = next((cc for k, cc in doc["cells"].items() if k.startswith("mixed_large|") and k.endswith(f"|{des}|all")), None)
                out.append(f"| {des} | {lab} | {t} | " + " | ".join(cls) + f" | {fmt(a[t]) if a else '–'} |")
    return "\n".join(out)


def block_confirm(runs):
    doc = runs.get("confirm_v3")
    if not doc:
        return "*Confirmatory run not available.*"
    c1, c2, c3 = [], [], []
    for key, c in doc["cells"].items():
        f, lv, des, label = key.split("|")
        for t in TESTS:
            d = c[t]
            if label == "all" and d["verdict"] == "EXCEEDS":
                c1.append(f"{f} {lv} {des} {t} ({fmt(d)})")
            if label == "all" and "fwer" in d and d["fwer"] * d["runs"] > 5:
                c2.append(f"{des} {t}: {round(d['fwer'] * d['runs'])} of {d['runs']} runs")
            if f == "mixed_large" and label != "all":
                lo, hi = (0.0, 0.015) if t == "shift" else (0.007, 0.014)
                if not (lo <= d["p01_unchanged"] <= hi):
                    c3.append(f"{des} {CLASS_NAMES[label]} {t}: {d['p01_unchanged']:.4f} (criterion [{lo}, {hi}])")
    res = lambda x: "**pass**" if not x else "**fail**: " + "; ".join(x)  # noqa: E731
    return "\n".join(["| criterion | result |", "|---|---|",
                      f"| C1 whole list: no EXCEEDS | {res(c1)} |",
                      f"| C2 complete null: ≤ 5 of 100 runs with any call | {res(c2)} |",
                      f"| C3 per-class tails (mixed_large): p01 in [0.007, 0.014] (shape_kernel, shape_any), ≤ 0.015 (shift) | {res(c3)} |"])


BLOCKS = dict(provenance=block_provenance, wholelist=block_wholelist, perclass=block_perclass,
              fdrstrata=block_fdrstrata, confirm=block_confirm)


def main(argv=None) -> int:
    argv = argv or sys.argv[1:]
    runs_dir, targets = Path(argv[0]), [Path(a) for a in argv[1:]]
    runs = {d.name: doc for d in sorted(runs_dir.iterdir()) if d.is_dir() and (doc := load(d))}
    for t in targets:
        text = t.read_text()
        for name, fn in BLOCKS.items():
            pat = re.compile(rf"(<!-- gen:{name} -->\n).*?(<!-- /gen:{name} -->)", re.S)
            text = pat.sub(lambda m, fn=fn: m.group(1) + fn(runs) + "\n" + m.group(2), text)
        t.write_text(text)
        print(t)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
