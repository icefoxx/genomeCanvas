"""Runtime and peak memory of ProfileSeg components, mapped onto the version timeline.

Each (regions, replicates) cell runs in a fresh process (single BLAS thread). Peak memory is the
tracemalloc peak of each component (numpy allocations are traced), plus the process RSS.

    v1 / v1.1   analyse (global; bootstrap 1000)
    v2          + per-region amplitude (moderated t, relative offsets) + moment shape tests
    v2.1        + shape_kernel + template_fit
    v2.2        + shift_score_test
    v2.3        + replicate_influence (leave-one-out, offsets)       [current default]
    options     shape_classes, align_profiles

    python benchmarks/profileseg_resources.py --out DIR
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
import tracemalloc
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

GRID = [(r, n) for r in (1000, 5000, 20000) for n in (2, 3, 4)]
KERNEL_ONLY = [(5000, 5), (5000, 6)]  # assignments grow fast: 5 vs 5 → 126, 6 vs 6 → 462


def measure(fn):
    tracemalloc.start()
    t = time.perf_counter()
    out = fn()
    dt = time.perf_counter() - t
    peak = tracemalloc.get_traced_memory()[1] / 2**20
    tracemalloc.stop()
    return out, dt, peak


def cell(regions: int, n: int, full: bool) -> dict:
    from genomecanvas.profileseg import metrics as M
    from genomecanvas.profileseg.analysis import analyse
    from genomecanvas.profileseg.data import Config
    from genomecanvas.profileseg.regional import (KERNEL_S0_QUANTILE, KERNEL_SIGMA_BP, TEMPLATE_SMOOTH,
                                                   amplitude_matrices, kernel_shape_test, moderated_t,
                                                   relative_offset, replicate_influence, shift_score_test,
                                                   template_fit)
    from genomecanvas.profileseg.shapeclass import align_profiles, shape_classes
    from genomecanvas.profileseg.simulate import SimSettings, simulate

    prof, _ = simulate("shift", SimSettings(n_regions=regions, n_ref=n, n_alt=n), seed=1)
    cfg = Config(reference="WT", alternative="KO")
    res = {}
    ref, alt = np.arange(n), np.arange(n, 2 * n)
    if full:
        _, res["global_bootstrap1000"], res["mem_global"] = measure(lambda: analyse(prof, cfg))
    a, res["global_no_bootstrap"], _ = measure(lambda: analyse(prof, cfg, bootstrap=0))

    def amplitude():
        comp, amp = amplitude_matrices(a)
        off, _ = relative_offset(amp["log2_signal"], ref, alt)
        amp["log2_signal_relative"] = amp["log2_signal"] - off
        return [moderated_t(amp[k], ref, alt) for k in amp]
    comp = np.flatnonzero(a.complete)
    sidx = comp[a.core_mass[:, comp].mean(0) > 0]
    y = a.core_signal[:, sidx]
    e = a.core_edges
    c = M.centers(e)
    core = M.core_mask(prof.edges, cfg.core_flank)
    cov = np.nan_to_num(a.values[:, sidx][..., core])
    if full:
        _, res["region_amplitude"], res["mem_region_amplitude"] = measure(amplitude)
        p_full = M.normalise(np.clip(y, 0, None) * M.widths(e))
        _, res["region_moments"], _ = measure(lambda: [moderated_t(M.centroid(p_full, e), ref, alt),
                                                      moderated_t(M.mean_abs_deviation(p_full, e), ref, alt)])
    _, res["shape_kernel"], res["mem_shape_kernel"] = measure(
        lambda: kernel_shape_test(y, np.ones(y.shape[1:], bool), ref, alt, KERNEL_SIGMA_BP, KERNEL_S0_QUANTILE, c))
    fit, res["template_fit"], res["mem_template_fit"] = measure(lambda: template_fit(y, cov, e, smooth=TEMPLATE_SMOOTH))
    _, res["shift_score"], _ = measure(lambda: shift_score_test(y, ref, alt, c, smooth=TEMPLATE_SMOOTH))
    if full:
        calls = {"shape_kernel": np.array(["Changed"] * 30 + ["Not_significant"] * (len(sidx) - 30)),
                 "shift": np.array(["Higher"] * 30 + ["Not_significant"] * (len(sidx) - 30))}
        _, res["replicate_influence"], res["mem_influence"] = measure(
            lambda: replicate_influence(y, cov, e, ref, alt, calls, 0.05, True))
        yp = np.clip(y.mean(0), 0, None); P = yp / np.maximum(yp.sum(-1, keepdims=True), 1e-12)
        _, res["shape_classes"], _ = measure(lambda: shape_classes(P, 8))
        _, res["align_profiles"], _ = measure(lambda: align_profiles(prof, cfg.background_width, cfg.core_flank, 150.0))
    res["tested_regions"] = int(len(sidx))
    res["max_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20  # bytes on macOS
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--cell", nargs=3, help=argparse.SUPPRESS)
    args = ap.parse_args(argv)
    if args.cell:
        r, n, full = int(args.cell[0]), int(args.cell[1]), args.cell[2] == "full"
        print(json.dumps(cell(r, n, full)))
        return 0
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for (r, n), full in [(g, True) for g in GRID] + [(g, False) for g in KERNEL_ONLY]:
        proc = subprocess.run([sys.executable, __file__, "--out", str(out), "--cell", str(r), str(n),
                               "full" if full else "kernel"], capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(proc.stderr[-2000:])
        d = json.loads(proc.stdout.strip().splitlines()[-1])
        d.update(regions=r, n=n)
        rows.append(d)
        print(f"done {r} regions {n}v{n}", flush=True)
    (out / "resources.json").write_text(json.dumps(rows, indent=1))
    lines = ["# ProfileSeg runtime and memory by component (simulation, 100 core bins, single thread)", "",
             "| regions | design | global (boot 1000) | region amplitude | kernel | template fit | shift score | influence | classes | align | peak RSS (MB) |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    f = lambda d, k: f"{d[k]:.2f}" if k in d else "—"
    for d in rows:
        lines.append(f"| {d['regions']} | {d['n']}v{d['n']} | {f(d, 'global_bootstrap1000')} | {f(d, 'region_amplitude')} | "
                     f"{f(d, 'shape_kernel')} | {f(d, 'template_fit')} | {f(d, 'shift_score')} | {f(d, 'replicate_influence')} | "
                     f"{f(d, 'shape_classes')} | {f(d, 'align_profiles')} | {d['max_rss_mb']:.0f} |")
    lines += ["", "## Cumulative time by version (s)", "",
              "| regions | design | v1/v1.1 | v2 | v2.1 | v2.2 | v2.3 (current default) |", "|---|---|---|---|---|---|---|"]
    for d in rows:
        if "global_bootstrap1000" not in d:
            continue
        v1 = d["global_bootstrap1000"]
        v2 = v1 + d["region_amplitude"] + d["region_moments"]
        v21 = v2 + d["shape_kernel"] + d["template_fit"]
        v22 = v21 + d["shift_score"]
        v23 = v22 + d["replicate_influence"]
        lines.append(f"| {d['regions']} | {d['n']}v{d['n']} | {v1:.1f} | {v2:.1f} | {v21:.1f} | {v22:.1f} | {v23:.1f} |")
    (out / "RESOURCES.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
