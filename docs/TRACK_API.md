# Track API

Track types understood by `render()`:

| `type` | Source | Role |
|---|---|---|
| `signal` | BedGraph, BigWig, or peak file | Filled step density |
| `overlay` | `members` track IDs | Line overlay of several signals |
| `peak` | BED | Interval ticks |
| `peak_signal` | narrowPeak/broadPeak | UCSC score-shaded blocks |
| `gene` | GTF, or figure `gene_model` | Transcript/exon models |
| `chromhmm` | 4-column BED | Chromatin-state bar |
| `cytoband` | UCSC cytoBandIdeo | Ideogram with view marker |
| `annotation` | BED | Generic intervals |

Region may be `chrom/start/end` or `gene` plus `gene_model` (or `query_ucsc: true`).

Python helpers that mirror trackplot's profile/analysis flow:

- `extract_profile` → `profile_summarize` → `profile_plot` / `profile_heatmap`
- `extract_summary` → `pca_scores` / `pca_plot`
- `diffpeak` → `volcano_plot`
- `summarize_signals` for condition collapse of locus slices

## BigWig engines

Switch backends with any of:

1. Function argument: `extract_profile(..., engine="bwtool")` or `query_signal(..., engine="pybigwig")`
2. Figure YAML: `bigwig_engine: bwtool`
3. CLI: `gcanvas render figure.yaml --bigwig-engine bwtool`
4. Environment: `GENOMECANVAS_BIGWIG_ENGINE=bwtool`
5. Process override: `set_bigwig_engine("pybigwig")` or `using_bigwig_engine("bwtool")`

Values: `auto`, `pybigwig`, `pybigtools`, `bwtool`.

`auto` uses **bwtool** for multi-region matrices when it is on PATH, and **pybigtools / pyBigWig** for single-locus queries (avoiding a subprocess per window). Explicit `engine=` always wins.

## Heatmap pipeline (deepTools-free)

`gcanvas heatmap` replaces `deeptools_heatmaps.wrapper.py`:

- region or genome-bin input
- reference-point (`center`/`TSS`/`TES`) or scale-regions (`--bodyLength`)
- per-`metaGrp` `profileMinMax` / `qnorm`
- Matplotlib heatmap + profile
- Pearson/Spearman, PCA, PairGrid
- kmeans/hclust/Leiden clustering
- voom differential BEDs
- peak-center Mann-Whitney between BED groups

See `examples/heatmap_pipeline/README.md`.

## Style, track heights, and future click edits

Figure-wide typography lives in `style:` (family, italic/bold, sizes, colors). Per-track
overrides use `tracks[].style` (`font_size`, `font_style`, `font_color`, `color`,
`fill_alpha`, `linewidth`, `height`). `tracks[].height` still sets relative panel height;
omitting it uses type defaults (peaks short, isoforms taller).

Drawn track types: `signal`, `derived_signal`, `difference`, `overlay`, `heatmap`,
`peak`, `peak_signal`, `gene`, `isoform`, `chromhmm`, `cytoband`, `annotation`,
`motif`, `loop`, `variant`.

Overlays (`overlays:`) are named labels/highlights/arrows. `render()` writes
`*.hits.json` pick targets. Click UIs should emit `FigureEdit` and call
`apply_edit(spec, edit)` then re-render — no matplotlib callbacks inside the engine.
