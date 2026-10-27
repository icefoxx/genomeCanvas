"""Label-free replicate structure versus the calls made in real-null splits (yeast insertion data).

For each genotype with four replicates, the distance between per-sample meta-profiles (mean over
regions of the background-subtracted, mass-normalised profile) is compared within and across the three
possible pairings (1,2|3,4), (1,3|2,4), (1,4|2,3): ratio = mean across-pair distance / mean within-pair
distance. A ratio well above 1 means the replicates fall into two systematically different pairs.
The measure uses no ProfileSeg test and no condition labels; the pairing is fixed by replicate numbering.

Joined with the per-split results of `profileseg_insertion.py` (insertion_runs.json): if the calls made
in a "null" split track this structure, the split is not a valid null (the two groups differ), rather
than the pooled null being liberal everywhere.

    python benchmarks/profileseg_replicate_structure.py YEAST_CACHE.npz insertion_runs.json OUT.md
"""

from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

import numpy as np

PAIRINGS = {"1+2": ((0, 1), (2, 3)), "1+3": ((0, 2), (1, 3)), "1+4": ((0, 3), (1, 2))}


def pairing_ratios(x: np.ndarray) -> dict:
    x = np.clip(x - np.median(x[..., :10], axis=-1, keepdims=True), 0, None)
    p = x / np.maximum(x.sum(-1, keepdims=True), 1e-12)
    meta = p.mean(1)
    dist = lambda a, b: float(np.linalg.norm(meta[a] - meta[b]))  # noqa: E731
    out = {}
    for name, (g1, g2) in PAIRINGS.items():
        within = (dist(*g1) + dist(*g2)) / 2
        across = np.mean([dist(a, b) for a in g1 for b in g2])
        out[name] = across / within
    return out


def main(argv=None) -> int:
    argv = argv or sys.argv[1:]
    d = np.load(argv[0], allow_pickle=True)
    v, ids, groups = np.nan_to_num(d["values"]), list(d["sample_ids"]), json.loads(str(d["groups"]))
    runs = [r for r in json.loads(Path(argv[1]).read_text()) if r["data"] == "yeast"]
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in runs:
        for m in ("shape_kernel", "shape_any", "shift"):
            key = "null" if r["insertion"] == "none" else ("amp" if r["insertion"] == "amplitude" else "fdp")
            by[r["split"]][f"{m}:{key}"].append(r["results"][m]["called" if key != "fdp" else "fdp"])
    lines = ["# Replicate structure versus calls in yeast real-null splits", "",
             "ratio = mean across-pair / mean within-pair distance of per-sample meta-profiles (label-free). "
             "Null calls: share of regions called with no insertion; FDP: mean over shift/broadening insertions.", "",
             "| genotype | split | pairing ratio | shape_kernel null calls | shape_any null calls | shift null calls | FDP kernel / any / shift |",
             "|---|---|---|---|---|---|---|"]
    rows = []
    for gt, smp in groups.items():
        ratios = pairing_ratios(v[[ids.index(s) for s in smp]])
        for name in PAIRINGS:
            split = f"{gt}:{smp[0]}+{smp[int(name[-1]) - 1]}"
            b = by.get(split)
            if not b:
                continue
            m = lambda k: float(np.mean(b[k]))  # noqa: E731
            rows.append((ratios[name], m("shape_kernel:null"), m("shape_kernel:fdp")))
            lines.append(f"| {gt} | {name.replace('+', '+')} vs rest | {ratios[name]:.2f} | {100 * m('shape_kernel:null'):.2f} % | "
                         f"{100 * m('shape_any:null'):.2f} % | {100 * m('shift:null'):.2f} % | "
                         f"{m('shape_kernel:fdp'):.3f} / {m('shape_any:fdp'):.3f} / {m('shift:fdp'):.3f} |")
    r = np.array(rows)
    hi = r[:, 0] > 1.2
    lines += ["", f"Splits with pairing ratio > 1.2: {int(hi.sum())} of {len(r)}; their mean `shape_kernel` FDP "
              f"{r[hi, 2].mean():.3f}, null calls in {int((r[hi, 1] > 0).sum())} of them. Other splits: mean FDP "
              f"{r[~hi, 2].mean():.3f}, null calls in {int((r[~hi, 1] > 0).sum())} of {int((~hi).sum())}. Spearman "
              f"(ratio, FDP) = {_spearman(r[:, 0], r[:, 2]):.2f}.", ""]
    Path(argv[2]).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


def _spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


if __name__ == "__main__":
    raise SystemExit(main())
