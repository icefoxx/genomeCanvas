---
name: genomecanvas-enrichment
description: GenomeCanvas region enrichment figures in the style of deepTools computeMatrix/plotHeatmap/plotProfile - signal heatmaps, mean profiles, profile + heatmap and k-means/hierarchical clustered heatmaps around region centres, TSS/TES or scaled region bodies.
---

# Enrichment & heatmaps

Run with the entry skill's runner ([genomecanvas-plot](../genomecanvas-plot/SKILL.md)).

| Plot slug | Draws |
|---|---|
| `enrichment/heatmap` | one heatmap panel per signal |
| `enrichment/profile` | mean signal per bin, one line per signal (± SEM) |
| `enrichment/profile-heatmap` | profile on top of each heatmap panel |
| `enrichment/clustered-heatmap` | heatmap with rows clustered (k-means by default) |

Inputs: ≥1 signal and ≥1 region file, at most 16 signals and 20,000 regions in total. Each region file
becomes a row group. Minus-strand regions (BED column 6) are flipped like in deepTools.

## Matrix settings (shared by all four)

| Key | Default | Notes |
|---|---|---|
| `alignment` | `center` | `center`, `TSS` (start), `TES` (end), `scale-regions` |
| `flank` | 1000 | bp either side, 20–1,000,000 |
| `binsize` | 50 | bp; `2 × flank / binsize` ≤ 4000 bins |
| `body_length` | 1000 | scaled body length, `scale-regions` only |
| `skip_zeros` | true | drop regions with no signal |

The runner applies the UI defaults (1000 / 50). The engine's own defaults (500 / 20) apply only if
you call the engine directly. Matrices are cached by input files and these settings, so restyling is
fast.

## Heatmap settings

| Key | Default | Notes |
|---|---|---|
| `sort` / `sort_by` | `descend` / `mean` | `ascend`, `none`; `median max min sum` |
| `cluster` / `clusters` | `none` (`kmeans` for clustered-heatmap) / 3 | `kmeans`, `hclust`; k 2–20, seed 7201 |
| `normalization` | `raw` | `profileMinMax` (per sample), `qnorm` (quantile) |
| `autoscale`, `zmin`, `zmax` | true | set `autoscale: false` to use zmin/zmax |
| `indi` | true | separate colour key per sample |
| `cmap` | `auto` | white → each sample's colour; or a Matplotlib name, or stops `#ffffff,#fdae61,#7f0000`; per-file `cmap` wins |
| `show_colorbar`, `show_group_labels` | true | |
| `interpolation` | `bilinear` | `nearest` for crisp pixels |
| `width_per_sample` | 2.7 in | `fig_width` overrides it |
| `colorbar_label`, `sample_label_size` | `signal`, 10 | |

## Profile settings

| Key | Default | Notes |
|---|---|---|
| `plot_scale` | `linear` | `log2` = log₂(x+1), `zscore` |
| `x_unit` | `auto` | `bp kb Mb` |
| `error_estimates` | true | SEM ribbon |
| `ln_v` | true | vertical guide at the anchor |
| `ln_h`, `ln_h_pos` | false, 0 | horizontal guide (`enrichment/profile` only) |
| `ylim_manual`, `ymin`, `ymax` | false | fixed y range |
| `legend`, `legend_pos` | true, `best` | `enrichment/profile` only |
| `xlab`, `ylab` | auto, `regions` | |
| `fill_alpha` | 0.18 | ribbon opacity |

## Outputs and pitfalls

* Writes `matrix.npz` and `matrix.json` next to the figure, so you can re-plot without recomputing.
* `cluster` reorders rows; the cluster labels C1…Ck replace region-file groups in the strip.
* `profileMinMax` and `qnorm` change what the colour scale means. Mention it whenever these figures
  are compared across samples.

Examples: `skills/genomecanvas-plot/examples/enrichment-*.json`.
