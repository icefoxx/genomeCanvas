# Changelog

## Unreleased

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
