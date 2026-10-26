"""Wall time and peak memory of ProfileSeg vs IWTomics, MMDiff2 and edgeR on the same simulated data.

One simulated experiment per size (3 vs 3, shift scenario). Each external tool runs alone in its
own Rscript process (timed with /usr/bin/time -l, single thread); ProfileSeg runs end to end
(analyse with 1000 bootstraps + per_region_tests) in a fresh Python process.

    python benchmarks/profileseg_tool_resources.py --out DIR [--sizes 300 1000]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import profileseg_external as EX
import profileseg_region_methods as RM
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.simulate import SimSettings

# IWTomics part of the comparison script only (MMDiff2 not loaded)
IWT_R = EX.R_SCRIPT.split("  # MMDiff2")[0].replace("library(MMDiff2); ", "") + "\n}\n"
MMD_R = r"""
suppressMessages({library(MMDiff2); library(GenomicRanges); library(Rsamtools)})
run <- commandArgs(TRUE)[1]
for (s in list.files(run, pattern = "\\.sam$", full.names = TRUE))
  asBam(s, sub("\\.sam$", "", s), overwrite = TRUE, indexDestination = TRUE)
meta <- list(ExpData = list(dataDir = run, sampleSheet = "samples.csv"))
peaks <- read.table(file.path(run, "peaks.tsv"), header = TRUE)
gr <- GRanges(peaks$chrom, IRanges(peaks$start + 1, peaks$end)); names(gr) <- peaks$name
md <- DBAmmd(meta); md <- setRegions(md, gr)
md <- suppressMessages(getPeakReads(md, pairedEnd = FALSE, run.parallel = FALSE))
md <- suppressMessages(estimateFragmentCenters(md))
md <- compHists(md, bin.length = 20, whichPos = "Center")
md <- suppressMessages(compDists(md, dist.method = "MMD", run.parallel = FALSE))
md <- setContrast(md, contrast = "byCondition")
md <- suppressMessages(compPvals(md, dist.method = "MMD"))
cat("ok\n")
"""


def timed(cmd: list[str]) -> dict:
    proc = subprocess.run(["/usr/bin/time", "-l"] + cmd, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(proc.stderr[-2000:])
    real = float(re.search(r"([\d.]+) real", proc.stderr).group(1))
    rss = int(re.search(r"(\d+)\s+maximum resident set size", proc.stderr).group(1)) / 2**20
    return dict(seconds=real, peak_rss_mb=rss)


PS_CODE = r"""
import os, sys
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ.setdefault(v, "1")
from genomecanvas.profileseg.simulate import SimSettings, simulate
from genomecanvas.profileseg.analysis import analyse
from genomecanvas.profileseg.data import Config
from genomecanvas.profileseg.regional import per_region_tests
prof, _ = simulate("shift", SimSettings(n_regions=int(sys.argv[1])), seed=int(sys.argv[2]))
a = analyse(prof, Config(reference="WT", alternative="KO"))
per_region_tests(a)
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--sizes", type=int, nargs="+", default=[300, 1000])
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="tools_", dir=out))
    (work / "iwt.R").write_text(IWT_R)
    (work / "mmd.R").write_text(MMD_R)
    (work / "edger.R").write_text(RM.R_SCRIPT)
    (work / "ps.py").write_text(PS_CODE)
    rows = []
    for size in args.sizes:
        seed = 4242
        run_dir = work / f"run_{size}"
        EX.prepare_run(run_dir, "shift", SimSettings(n_regions=size), Config(reference="WT", alternative="KO"), seed)
        d = RM.prepare("shift", SimSettings(n_regions=size), None, seed)
        RM.count_tables(d, RM.support(d), run_dir)
        row = dict(regions=size)
        row["ProfileSeg (global + per region)"] = timed([sys.executable, str(work / "ps.py"), str(size), str(seed)])
        row["IWTomics (4 tests, B = 1000)"] = timed(["Rscript", str(work / "iwt.R"), str(run_dir)])
        row["edgeR QL diffSplice (2 tables)"] = timed(["Rscript", str(work / "edger.R"), str(run_dir)])
        row["MMDiff2 (BAM → p-values)"] = timed(["Rscript", str(work / "mmd.R"), str(run_dir)])
        rows.append(row)
        print(json.dumps(row), flush=True)
    (out / "tool_resources.json").write_text(json.dumps(rows, indent=1))
    tools = [k for k in rows[0] if k != "regions"]
    lines = ["# Wall time and peak memory, 3 vs 3, single thread", "",
             "| tool | " + " | ".join(f"{r['regions']} regions: time (s) / peak RSS (MB)" for r in rows) + " |",
             "|---|" + "---|" * len(rows)]
    for t in tools:
        lines.append(f"| {t} | " + " | ".join(f"{r[t]['seconds']:.1f} / {r[t]['peak_rss_mb']:.0f}" for r in rows) + " |")
    (out / "TOOL_RESOURCES.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
