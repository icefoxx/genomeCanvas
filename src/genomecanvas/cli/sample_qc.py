"""gcanvas-sampleqc: sample QC for large ChIP-seq collections.

    gcanvas-sampleqc matrix  --sheet samples.tsv --sample-col sample_id --path-col bw_path_hpc \\
                             --out matrix.npz [--bin 10000] [--workers 16] [--blacklist bl.bed]
    gcanvas-sampleqc preview --sheet samples.tsv [--matrix matrix.npz]
    gcanvas-sampleqc run     --sheet samples.tsv --matrix matrix.npz --out qc/ \\
                             [--mark target --group model --batch batch] [--filter "status=PASS"]

`matrix` is the slow step (reads every BigWig once; resumable). `run` takes seconds and can be
repeated with other columns, filters and thresholds; Studio's "Sample QC" plot runs the same code.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _options(args) -> dict:
    return {
        "qc_sample": args.sample_col,
        "qc_mark": args.mark,
        "qc_group": args.group,
        "qc_batch": args.batch,
        "qc_filters": "; ".join(args.filter or []),
        "qc_k": args.k,
        "qc_top": args.top,
        "qc_method": args.method,
        "qc_z": args.z,
        "qc_margin": args.margin,
        "qc_skip_marks": args.skip_marks,
    }


def cmd_matrix(args) -> int:
    from ..analysis.sample_qc import apply_filters, parse_filters, read_sheet
    from ..pipeline.sample_matrix import build_matrix

    columns, rows = read_sheet(args.sheet)
    for column in (args.sample_col, args.path_col):
        if column not in columns:
            raise SystemExit(
                f"Column {column!r} is not in {args.sheet}; columns: {', '.join(columns)}"
            )
    rows = apply_filters(rows, parse_filters("; ".join(args.filter or [])))
    sources = {}
    for row in rows:
        path = row[args.path_col]
        for mapping in args.path_map or []:
            old, _, new = mapping.partition("=")
            if path.startswith(old):
                path = new + path[len(old) :]
        if row[args.sample_col] in sources:
            raise SystemExit(f"Duplicate sample {row[args.sample_col]!r} in {args.sample_col}")
        sources[row[args.sample_col]] = path
    chroms = [c.strip() for c in (args.chroms or "").split(",") if c.strip()] or None
    report = build_matrix(
        sources,
        args.out,
        size=args.bin,
        stat=args.stat,
        chroms=chroms,
        blacklist=args.blacklist,
        workers=args.workers,
        checkpoint=args.checkpoint,
        log=lambda line: print(line, file=sys.stderr, flush=True),
    )
    print(json.dumps(report, indent=2))
    return 1 if report["missing_files"] or report["failed"] else 0


def cmd_preview(args) -> int:
    from ..analysis.sample_qc import preview

    out = preview(args.sheet, args.matrix, _options(args))
    out.pop("rows")
    print(json.dumps(out, indent=2, default=str))
    return 0


def cmd_run(args) -> int:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    from ..analysis import sample_qc_plots as plots
    from ..analysis.sample_qc import qc_from_files, sorted_rows, write_tables

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    result, notes = qc_from_files(args.sheet, args.matrix, _options(args))
    write_tables(result, out)
    figure = Figure(figsize=(17, 5.6), layout="constrained")
    plots.draw_overview(result, figure)
    for ext in ("png", "pdf"):
        figure.savefig(out / f"qc_overview.{ext}", dpi=200)
    plots.pca_facets(result, out / "qc_pca_by_mark")
    plots.mark_pages(result, out / "qc_marks.pdf")
    counts: dict = {}
    for row in result.rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    flagged = [r for r in sorted_rows(result) if r["status"] not in ("ok", "undetermined")]
    print(
        json.dumps(
            dict(
                out=str(out),
                status_counts=counts,
                notes=notes,
                top_flags=[
                    {k: r.get(k) for k in ("sample", "mark", "group", "status", "reasons")}
                    for r in flagged[:20]
                ],
            ),
            indent=2,
            default=str,
        )
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gcanvas-sampleqc",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)
    matrix = sub.add_parser(
        "matrix", help="bin every sample's BigWig into one matrix (slow, resumable)"
    )
    matrix.add_argument("--path-col", required=True, help="sheet column with BigWig paths")
    matrix.add_argument(
        "--path-map", action="append", metavar="OLD=NEW", help="rewrite a path prefix"
    )
    matrix.add_argument("--out", required=True, help="output .npz")
    matrix.add_argument("--bin", type=int, default=10_000)
    matrix.add_argument("--stat", default="mean", choices=["mean", "mean0", "sum", "max"])
    matrix.add_argument(
        "--chroms", help="comma-separated chromosomes (default: all but chrM, random, …)"
    )
    matrix.add_argument("--blacklist", help="BED of regions to drop (e.g. ENCODE blacklist)")
    matrix.add_argument("--workers", type=int, default=4)
    matrix.add_argument(
        "--checkpoint", help="folder for per-sample checkpoints (default OUT.parts)"
    )
    preview = sub.add_parser(
        "preview", help="columns, guesses, matches and tiers of a sample sheet"
    )
    preview.add_argument("--matrix")
    run = sub.add_parser("run", help="flag table, overview, per-mark PCA and heatmaps")
    run.add_argument("--matrix", required=True)
    run.add_argument("--out", required=True)
    for p in (matrix, preview, run):
        p.add_argument("--sheet", required=True)
        p.add_argument(
            "--sample-col", default=None, help="sheet column matching matrix sample names"
        )
        p.add_argument(
            "--filter", action="append", help='e.g. "status=PASS" or "target!=INPUT,IgG"'
        )
    for p in (preview, run):
        p.add_argument("--mark")
        p.add_argument("--group")
        p.add_argument("--batch")
        p.add_argument("--k", type=int, default=5, help="neighbours voting on the mark")
        p.add_argument("--top", type=int, default=5000, help="most variable bins used")
        p.add_argument("--method", default="pearson", choices=["pearson", "spearman"])
        p.add_argument("--z", type=float, default=3.0, help="robust SDs for the low-quality floor")
        p.add_argument("--margin", type=float, default=0.02, help="r gap needed to call a mislabel")
        p.add_argument("--skip-marks", default="INPUT, IgG", help="marks without a group check")
    args = parser.parse_args(argv)
    if args.command == "matrix" and not args.sample_col:
        parser.error("matrix needs --sample-col")
    try:
        return {"matrix": cmd_matrix, "preview": cmd_preview, "run": cmd_run}[args.command](args)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
