---
name: genomecanvas-differential
description: GenomeCanvas exploratory differential figures - differential-peak volcano between two sample groups over BED regions, and peak-center comparison heatmaps between region sets. Use for a quick look; formal differential binding needs DESeq2/edgeR/DiffBind.
---

# Differential analysis (exploratory)

Run with the entry skill's runner ([genomecanvas-plot](../genomecanvas-plot/SKILL.md)).

## `differential/volcano`

Inputs: ≥4 signals, ≥1 region file. Each signal's `group` defines its condition. The numerator and
denominator groups need **≥2 replicates each**.

| Key | Default | Notes |
|---|---|---|
| `numerator` | `Treatment` | must equal a `group` value exactly |
| `denominator` | `Control` | different from `numerator` |
| `fdr` | 0.05 | threshold for colouring |
| `up_color`, `down_color` | palette | |
| `rel_stat`, `rel_skip_zeros`, `rel_max_rows` | `sum`, true, all | per-region summary, as in [relationships](../genomecanvas-relationships/SKILL.md) |

Regions are always the BED intervals; genome bins are not available here. Outputs `statistics.tsv`
(log2 fold change, p, FDR) and `summary.tsv`. The engine always adds a warning that the test uses a normal-tail
approximation. Pass that warning on, and do not present the result as a formal differential analysis.

## `differential/peak-center`

Inputs: ≥1 signal and ≥2 region sets. It compares signal at the centre of each region set
(Mann–Whitney, **unadjusted** p-values).

| Key | Default | Notes |
|---|---|---|
| `center_bins` | 5 | half-window in bins around the centre, 1–200 |
| `pc_cmap` | `viridis` | p-value colour scale |
| `alignment`, `flank`, `binsize`, `skip_zeros` | `center`, 1000, 50, true | matrix settings from [enrichment](../genomecanvas-enrichment/SKILL.md) |

Extra outputs: `comparison.peakCenter.pvalues.tsv` and `comparison.peakCenter.signals.clustermap.*`.

## Pitfalls

* A group-name mismatch gives "Volcano plots need two distinct groups with ≥2 replicates each" with
  the current group counts. Fix the `group` fields or `numerator`/`denominator`.
* Signal tracks are usually normalised coverage, not counts, so the fold changes are descriptive only.

Examples: `skills/genomecanvas-plot/examples/differential-*.json`.
