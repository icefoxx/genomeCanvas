# Real-data reference comparison

The inputs in `tests/fixtures/reference_demo/` are the public demo
BigWig/GRCm38 data. They are not stored in git; fetch them first with
`python tests/fixtures/reference_demo/download.py`. The comparison script uses only the public data files and
GenomeCanvas-side rendering code; it does not import or copy the reference
plotting implementation.

```bash
python examples/reference_comparison/generate.py
```

Outputs are written to this directory as `p01_genomecanvas_realdata.png` and
`.pdf`. They use the same four locus coordinates and four signal categories
visible in the reference faceted demo, while retaining independent layout and
rendering decisions. The script exposes `Y_SCALE_MODE`, `Y_LIMITS`,
`SHOW_GENE_TRACK`, `SHOW_CYTOBAND`, and `X_LABEL` near the top for direct
customization. Set a row to `(0, 400)` or `(0, 200)` for a fixed y-axis.
The p03 region renderer also accepts a `track_colors` list for independent
signal colors.
