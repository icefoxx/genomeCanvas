# Real BigWig / BED sources used by GenomeCanvas capability checks

Downloaded or copied from public repositories. Not synthetic mini-genome fixtures.

## Official deepTools Galaxy wrapper tests

Source: https://github.com/deeptools/deepTools/tree/master/galaxy/wrapper/test-data

- `official/computeMatrix2.bw` + `official/computeMatrix2.bed` — `computeMatrix` / `plotHeatmap` inputs
- `official/computeMatrix1.bed` — phiX174 intervals (pairs with bamCoverage outputs)
- `official/bamCoverage_result4.bw`, `official/bamCoverage_result5.bw` — `bamCoverage` example BigWigs (`phiX174`)

## pyGenomeTracks (deepTools org) Drosophila chrX

Source: https://github.com/deeptools/pyGenomeTracks/tree/master/pygenometracks/tests/test_data

- `bw/pgt_chrx_2e6_5e6.bw` ← `bigwig_chrx_2e6_5e6.bw` (chrom `chrX`)
- `bw/pgt_chrx_2.5e6_3.5e6.bw` ← `bigwig2_X_2.5e6_3.5e6.bw` (chrom `X`)
- `bw/pgt_chrx_2.5e6_3.5e6_b.bw` ← `bigwig3_X_2.5e6_3.5e6.bw` (chrom `X`)
- `bed/HoxD_cluster_regulatory_regions_mm10.bed`, `bed/HoxD.gtf`

## GEO mouse 2-cell (`tests/fixtures/reference_demo`)

GEO: GSM5746910–GSM5746913 (mm10 / GRCm38). Not stored in git; fetch with
`python tests/fixtures/reference_demo/download.py`.

## H1 ESC trackplot demo (human chr6 OCT4)

`tests/fixtures/trackplot_demo/` — PoisonAlien/trackplot demo BigWigs and peaks.

## Fly genes for scale-regions

`bed/dm3_chrX_2_5Mb_refGene.bed` — UCSC dm3 `refGene` chrX:2–5 Mb
(`https://api.genome.ucsc.edu/getData/track?genome=dm3;track=refGene;chrom=chrX;start=2000000;end=5000000`).
