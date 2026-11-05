---
name: genomecanvas-tracks
description: GenomeCanvas genome-browser style locus figures - stacked signal tracks, overlaid signals, and gene/transcript models with peaks, chromHMM states and highlight windows. Use when the user wants to show signal at one genomic locus.
---

# Genomic tracks

Run with the entry skill's runner ([genomecanvas-plot](../genomecanvas-plot/SKILL.md)).

| Plot slug | Inputs | Draws |
|---|---|---|
| `tracks/stacked` | ≥1 signal; optional BED, chromHMM BED, GTF | one panel per signal, then region tracks, then genes |
| `tracks/overlay` | ≥2 signals; optional BED, GTF | the signal panels plus one panel with all signals overlaid |
| `tracks/genes` | 1 GTF; optional BED | gene models only (signals are ignored) |

## Inputs

* Signals: BigWig or bedGraph, at most 16 per figure.
* Region files become peak tracks coloured like their file. A region file is drawn as a **chromHMM
  state track** when its name contains `chromhmm`, `hmm` or `segway`, or when every feature name
  looks like `E1`, `E12_TxWk`, ….
* Only the first GTF is used. If a GTF is present, region files whose name or label contains "gene"
  or "transcript" are skipped (to avoid drawing genes twice).
* Use `gcplot.py inspect FILE` to see the chromosomes and data extent before choosing a locus.

## Settings

| Key | Default | Notes |
|---|---|---|
| `chrom` | `chr1` | must match the files' naming |
| `start` / `end` | 0 / 100000 | 0-based, end-exclusive; window ≤ 20 Mb |
| `assembly` | | label only, e.g. `hg38` |
| `highlights` | | one per line: `start-end [colour] [opacity]`, e.g. `31136000-31138000 #f5c542 0.25`; commas and `chr6:` prefixes are accepted |
| `scale` | `auto_highest` | shared y range; `per_track` for independent ranges |
| `transcripts` | `collapsed` | `individual` shows every transcript |
| `gene_label` | `name` | `id` or `none` |
| `x_tick_unit` | `auto` | `bp kb Mb` |
| `track_gap` | 0.18 | 0–0.8 |
| `y_ticks`, `y_decimals` | 3, 0 | |
| `linewidth`, `fill_alpha` | 0.6, 0.85 | |
| `track_label_size` | 9 | |
| `xlab` | coordinate label | |

`advanced.layout`, `advanced.style` and `advanced.track_style` (applied to every signal track's
`style`) pass straight into the figure spec (`src/genomecanvas/model/figure_spec.py`) for anything
not listed here.

## Pitfalls

* Shared y scale (`auto_highest`) lets a strong track flatten weak ones. Choose it when tracks should be
  comparable and `per_track` when you want to show the shape of each track.
* An empty panel usually means the chromosome name is wrong (`chr6` vs `6`) or the window is outside
  the data. Run `inspect` to check.

Example: `skills/genomecanvas-plot/examples/tracks-stacked.json`.
