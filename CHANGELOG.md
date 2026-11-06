# Changelog

## Unreleased

- Heatmap colour (z) scale and profile y axis: new `auto` default groups panels whose ranges are within a set fold (`--zScale/--zScaleRatio`, `--profileYAxis/--profileShareRatio`; modes auto, shared, independent, group); `--zMin/--zMax` and new `--profileYMin/--profileYMax` accept one value or one per sample.
- Locus tracks: new `auto` y-scale mode (Studio default) shares an axis only between tracks within 3-fold; tick labels get enough decimals for small signals.
- Studio: gene labels now show gene names (the gene track ignored the label setting and showed IDs); manifest rows may be symlinks inside the project folder.
- Added Sample QC for large ChIP-seq collections (`gcanvas-sampleqc`, Studio "Sample relationships → Sample QC"): parallel, resumable bins × samples matrix builder; kNN mark-label check; per-mark group check with low-quality, mislabel, swap and batch-confounded flags; flag table, per-mark PCA and heatmap pages; Studio sample-sheet preview with column mapping and filters.
- Studio accepts table files (`.tsv`, `.csv`, `.npz`) as a new `table` role.
- Added `skills/`: agent skills for every Studio plot with one runner (`gcplot.py`) for local and server use.
- Added trackplot-parity locus tracks: overlay, peak, peak_signal, chromHMM, gene, cytoband, gene-symbol regions, and highlights.
- Added profile extract/summarize/plot/heatmap, condition collapse, PCA, differential peaks, and volcano plots.
- Added figure-wide and per-track font/color/italic styling, type default track heights, and overlay/hit/FigureEdit APIs for later click annotations.
- Wired isoform, loop, variant, motif, difference, and member heatmap tracks in the locus renderer.
- Added a deepTools-free heatmap pipeline (`gcanvas heatmap`) covering grouped matrices, renormalization, clustering, voom BEDs, bin mode, and scale-regions.
- Added selectable BigWig engines (`auto`, `pybigwig`, `pybigtools`, `bwtool`) via function argument, YAML, CLI, environment variable, and `using_bigwig_engine`.
- Added BigWig backends for pybigtools, pyBigWig, and optional bwtool.
- Added optional UCSC REST lookup for gene models and cytobands.
- Added the initial project scaffold (TASK-001).
- Added typed FigureSpec validation and YAML/JSON serialization.
- Added coordinate-safe readers, deterministic signal transforms, tracks, rendering, export, and provenance.
- Added TASK-031 local background-adjusted visualization transforms with explicit provenance and validation.
