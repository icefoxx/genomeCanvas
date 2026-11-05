---
name: genomecanvas-transform
description: GenomeCanvas signal transformations at one locus - difference, ratio, log2 ratio between two signals, and local background correction - drawn as A, B and the derived track. Use to compare two tracks (e.g. treatment vs control, ChIP vs input) in a region.
---

# Signal transformations

Run with the entry skill's runner ([genomecanvas-plot](../genomecanvas-plot/SKILL.md)).

Inputs: ≥2 signals. Only the **first two** in `files` are used: A = first, B = second. Set
`direction: "a-b"` to swap them. The figure has three panels: A, B and the derived track.
`signals.tsv` holds `start, A, B, derived` per bin.

| Plot slug | Derived track |
|---|---|
| `transform/difference` | A − B |
| `transform/ratio` | A / B (with pseudocount; dashed baseline at 1) |
| `transform/log2-ratio` | log₂((A + p) / (B + p)) |
| `transform/background` | A − baseline, where baseline = median (or mean) of **B** in the ±flank windows outside the locus |

## Settings

| Key | Default | Notes |
|---|---|---|
| `chrom`, `start`, `end` | `chr1`, 0, 100000 | locus; 0-based, end-exclusive |
| `pixels` | 800 | bins across the locus, 50–5000 |
| `pseudocount` | 0.01 | must be > 0 (ratio, log2 ratio) |
| `bg_flank` | 5000 | background only; the runner sends it as `flank`, as the UI does |
| `estimator` | `median` | `mean` (background only) |
| `direction` | `b-a` | `a-b` swaps A and B |
| `derived_color` | `#3b4252` | |
| `linewidth`, `fill_alpha` | 1.5, 0.25 | |

## Pitfalls

* Background correction takes the baseline from signal **B**'s flanks and subtracts it from A. With
  an input/IgG track as B this gives "A minus local input level". To use A's own flanks, list A
  twice.
* These are visualisation transforms (see `docs/LOCAL_BACKGROUND.md`). Call the result
  "background-adjusted signal" or "log₂ ratio", not "normalised".
* A small pseudocount makes ratios explode where B ≈ 0. Raise it, or prefer the difference.

Examples: `skills/genomecanvas-plot/examples/transform-*.json`.
