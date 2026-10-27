# ProfileSeg v3.1 — review packet

*Version 3.1.0, branch `profileseg-v3`. Consolidated after six review rounds; the full record (withdrawn claims,
superseded analyses, responses per round) is in `docs/PROFILESEG_V3_REVIEW_HISTORY.md`. Every table here is
generated from result files under `benchmarks/results/` (generation script: `benchmarks/profileseg_review_tables.py`
for the calibration blocks; the others were copied programmatically from the JSON files named in each section).
"Seg" means region-level screening of differential profiles, not segmentation of curves into intervals.*

## 1. Summary

ProfileSeg compares signal profiles around predefined anchors (TSS, summits, nucleosome dyads) between two
conditions with biological replicates, and reports amplitude, shape and a translation-type effect separately.
The questions it addresses: did the spatial distribution of signal at a region change, and, where a translation
describes the change, by how many bp. The statistical core is a per-region permutation null pooled across
regions, needed because 2 or 3 replicates per group give only 3 or 10 label assignments per region.

What the evidence supports, and what it does not:

| claim | status |
|---|---|
| whole-list FDR for the shape tests near or below 0.05 across peak width, depth, noise, share changed (0–50 %), overlapping windows, balanced artefacts and real replicate noise, 2 and 3 replicates | supported empirically; pre-specified confirmatory run on new seeds passed; not a proven bound |
| per-class tail calibration (narrow / broad / low-depth regions) within the pre-specified tolerance (≤ ~1.5× nominal) | supported with ≥ 3000 regions at 2 vs 2 and ≥ 1000 at 3 vs 3 in the tested mix; **not** with fewer (2 vs 2: up to 2.6× liberal for narrow peaks) |
| FDR within a subgroup (e.g. broad or shallow regions) under global BH | **not supported**; optional stratified BH targets it within matching strata and loses the combined-list level |
| exact permutation validity across regions | not claimed (the null is pooled; cross-region comparability is empirical) |
| real-data calibration | PRO-seq 3 vs 3 within nominal; **yeast MNase 2 vs 2 above nominal** (`shape_kernel` FDP 0.11–0.13, `shape_any` 0.10–0.11); the excess is associated with replicate structure, cause undetermined |
| more sensitive than edgeR diffSpliceDGE | pre-specified comparison met for 50 bp (narrow) and 100 bp (medium) translations at BH 5 %; at comparable observed FDP the advantage is robust for 50 bp on narrow peaks, small or absent for 100 bp on medium peaks, modest for broadening at 3 vs 3; neither detects broadening at 2 vs 2; simulation only |
| translation estimate (bp) | a consistency-checked shift that best matches the two condition-mean profiles; direction correct in every evaluated cell; estimate ÷ true 0.94–1.08 where `shift` power ≥ 0.3; biased upwards (up to 1.38) near the detection limit (selection) |
| ProfileSeg identifies translations as a mechanism | **not supported**: `translation_consistent` is a fit-consistency flag; one-sided gains that reach significance pass it; a stricter cut was evaluated and rejected |
| per-region absolute amplitude | **miscalibrated** under complete nulls in two designs (FWER = FDR 0.20 and 0.275); not claimed |

Positioning: a replicate-aware, region-level screening tool for shape changes with an empirically studied
calibration and explicit limits; not yet a validated tool for quantifying translations, and not shown to be
generally superior to count-based alternatives.

## 2. Method

- **Inputs.** Sample manifest (condition, replicate, optional batch) with bigWig or BAM signal files (BAM:
  fragment midpoints, coverage or read ends; MAPQ, duplicate, proper-pair and fragment-size filters), or a
  precomputed matrix; a BED of predefined anchors with strand; window, bin and background settings.
- **Shape tests** on background-corrected, mass-normalised profiles:
  `shape_kernel` (MMD-type kernel distance between condition means, scaled by within-condition spread plus a
  label-invariant s0); `shift` (projection of each sample on the pooled template's derivative; directional);
  `shape_any` (per-region Cauchy combination of the two p-values, then one BH — approximately valid through a
  tail approximation, calibration checked by simulation).
- **Null.** Statistics of all non-observed label assignments, pooled over regions; in v3 pooled within quantile
  strata of each region's total spread (identical under every assignment) when each stratum keeps ≥ 200
  regions and ≥ 2000 null values. The output reports the number of testable regions, strata and null values.
- **Translation description.** `estimate` = cross-correlation shift of the condition means (all lags). A
  direction (Higher/Lower) and bp value are given only when it agrees in sign with the projection, explains
  ≥ 50 % of the between-condition difference (`translation_r2`) and is unambiguous (`translation_ambiguity`
  < 0.9); otherwise `Changed`, NaN bp and the reason in `translation_flag`. `translation_fit` (high_r2 ≥ 0.8 /
  moderate_r2) is descriptive. The earlier per-sample template fit is kept as `shift_bp_fit`.
- **Diagnostics.** Depth per region (exact for BAM, estimated for bigWig), replicate influence, batch
  confounding (tests refused when condition and batch coincide), fragment-length QC for BAM.

## 3. Calibration of the shape tests

Script `benchmarks/profileseg_calibration_boundary.py`; details and deviations from protocol in
`docs/PROFILESEG_CALIBRATION.md`. Development runs reuse the seeds on which the v3 null was designed; the
confirmatory run used new seeds and criteria fixed beforehand.

### 3.1 Runs

<!-- gen:provenance -->
| run | purpose | seeds | null | FDR | commit | uncommitted src files | regional.py sha256 | src diff sha256 | started → finished |
|---|---|---|---|---|---|---|---|---|---|
| confirm_v3 | confirmatory, v3 default (frozen) | 9100000+ | strata | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T21:05:10 → 2026-10-26T21:17:20 |
| dev_fdr_depth | development, depth FDR strata | 5000000+ | strata | depth | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T21:01:48 → 2026-10-26T21:05:09 |
| dev_fdr_spread | development, spread FDR strata | 5000000+ | strata | spread | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:58:23 → 2026-10-26T21:01:46 |
| dev_v24 | development, v2.4 null | 5000000+ | none | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:46:07 → 2026-10-26T20:58:22 |
| dev_v3 | development, v3 default | 5000000+ | strata | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:33:37 → 2026-10-26T20:46:06 |
<!-- /gen:provenance -->

The confirmatory run was repeated with the current translation code: all 50 original cells identical
(`benchmarks/results/round4/confirm_v3_rerun`).

### 3.2 Whole-list FDR and confirmatory criteria

<!-- gen:wholelist -->
| run | cells × tests judged | ok | WATCH | EXCEEDS | largest mean FDP | largest CI upper bound |
|---|---|---|---|---|---|---|
| dev_v3 | 114 | 114 | 0 | 0 | 0.052 | 0.059 |
| dev_v24 | 114 | 114 | 0 | 0 | 0.052 | 0.063 |
| confirm_v3 | 114 | 114 | 0 | 0 | 0.049 | 0.058 |

Selected cells (FDP [95 % CI] / power; whole list):

| factor | level | design | run | shape_kernel | shift | shape_any |
|---|---|---|---|---|---|---|
| width | 60 | 2v2 | dev_v3 | 0.052 [0.046, 0.057] / 0.87 | 0.034 [0.029, 0.040] / 0.48 | 0.046 [0.041, 0.050] / 0.83 |
| width | 60 | 2v2 | confirm_v3 | 0.046 [0.041, 0.051] / 0.86 | 0.034 [0.028, 0.039] / 0.48 | 0.039 [0.035, 0.044] / 0.81 |
| width | 400 | 3v3 | dev_v3 | 0.027 [0.015, 0.039] / 0.02 | 0.011 [0.000, 0.024] / 0.00 | 0.014 [0.004, 0.024] / 0.01 |
| width | 400 | 3v3 | confirm_v3 | 0.032 [0.017, 0.047] / 0.03 | 0.000 [0.000, 0.000] / 0.00 | 0.006 [0.000, 0.012] / 0.01 |
| depth | 5.0 | 3v3 | dev_v3 | 0.037 [0.031, 0.044] / 0.26 | 0.036 [0.025, 0.047] / 0.07 | 0.030 [0.023, 0.037] / 0.22 |
| depth | 5.0 | 3v3 | confirm_v3 | 0.045 [0.037, 0.052] / 0.28 | 0.043 [0.029, 0.057] / 0.07 | 0.041 [0.034, 0.048] / 0.24 |
| noise | high | 3v3 | dev_v3 | 0.036 [0.028, 0.043] / 0.18 | 0.029 [0.014, 0.044] / 0.01 | 0.034 [0.025, 0.042] / 0.14 |
| noise | high | 3v3 | confirm_v3 | 0.034 [0.026, 0.042] / 0.19 | 0.026 [0.011, 0.040] / 0.02 | 0.032 [0.024, 0.041] / 0.13 |
| proportion | 0.0 | 2v2 | dev_v3 | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] |
| proportion | 0.0 | 2v2 | confirm_v3 | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] |
| proportion | 0.0 | 3v3 | dev_v3 | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] |
| proportion | 0.0 | 3v3 | confirm_v3 | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] | FWER 0.000 [0.000, 0.036] |
| proportion | 0.5 | 2v2 | dev_v3 | 0.043 [0.040, 0.046] / 0.63 | 0.030 [0.026, 0.033] / 0.32 | 0.037 [0.034, 0.040] / 0.59 |
| proportion | 0.5 | 2v2 | confirm_v3 | 0.044 [0.040, 0.047] / 0.63 | 0.031 [0.027, 0.035] / 0.31 | 0.039 [0.036, 0.043] / 0.59 |
| correlation | overlap | 2v2 | dev_v3 | 0.046 [0.041, 0.051] / 0.45 | 0.043 [0.035, 0.051] / 0.23 | 0.045 [0.039, 0.051] / 0.43 |
| correlation | overlap | 2v2 | confirm_v3 | 0.048 [0.042, 0.054] / 0.44 | 0.045 [0.037, 0.053] / 0.22 | 0.047 [0.041, 0.053] / 0.42 |
| batch | balanced 15 bp | 3v3 | dev_v3 | 0.031 [0.027, 0.035] / 0.62 | 0.029 [0.024, 0.034] / 0.36 | 0.030 [0.026, 0.034] / 0.57 |
| batch | balanced 15 bp | 3v3 | confirm_v3 | 0.030 [0.026, 0.034] / 0.61 | 0.027 [0.021, 0.032] / 0.35 | 0.030 [0.025, 0.034] / 0.57 |
| semireal | real noise | 2v2 | dev_v3 | 0.049 [0.043, 0.055] / 0.54 | 0.038 [0.031, 0.045] / 0.42 | 0.048 [0.042, 0.054] / 0.56 |
| semireal | real noise | 2v2 | confirm_v3 | 0.047 [0.039, 0.055] / 0.53 | 0.042 [0.033, 0.051] / 0.42 | 0.048 [0.040, 0.057] / 0.55 |
| mixed_large | 3 classes × 1000 | 2v2 | dev_v3 | 0.045 [0.042, 0.049] / 0.32 | 0.037 [0.032, 0.042] / 0.14 | 0.040 [0.036, 0.043] / 0.29 |
| mixed_large | 3 classes × 1000 | 2v2 | confirm_v3 | 0.044 [0.041, 0.048] / 0.32 | 0.034 [0.029, 0.038] / 0.13 | 0.037 [0.034, 0.040] / 0.28 |
| mixed_large | 3 classes × 1000 | 3v3 | dev_v3 | 0.045 [0.042, 0.047] / 0.44 | 0.042 [0.038, 0.045] / 0.19 | 0.043 [0.040, 0.046] / 0.41 |
| mixed_large | 3 classes × 1000 | 3v3 | confirm_v3 | 0.046 [0.043, 0.048] / 0.44 | 0.039 [0.036, 0.043] / 0.19 | 0.043 [0.040, 0.046] / 0.41 |
<!-- /gen:wholelist -->

<!-- gen:confirm -->
| criterion | result |
|---|---|
| C1 whole list: no EXCEEDS | **pass** |
| C2 complete null: ≤ 5 of 100 runs with any call | **pass** |
| C3 per-class tails (mixed_large): p01 in [0.007, 0.014] (shape_kernel, shape_any), ≤ 0.015 (shift) | **pass** |
<!-- /gen:confirm -->

"ok" means no statistically clear exceedance (the 95 % Monte Carlo CI of the mean FDP not above 0.05), not a
demonstrated bound.

### 3.3 Per-class tails

<!-- gen:perclass -->
Share of **unchanged** regions with p < 0.01 (nominal 0.01), mixed designs:

| design | class | test | v2.4 null (dev) | v3 null (dev) | v3 null (confirmatory) |
|---|---|---|---|---|---|
| mixed 2v2 | narrow (sd 60 bp) | shape_kernel | 0.0225 | 0.0225 | 0.0257 |
| mixed 2v2 | narrow (sd 60 bp) | shift | 0.0023 | 0.0023 | 0.0029 |
| mixed 2v2 | narrow (sd 60 bp) | shape_any | 0.0131 | 0.0131 | 0.0136 |
| mixed 2v2 | broad (sd 400 bp) | shape_kernel | 0.0000 | 0.0000 | 0.0000 |
| mixed 2v2 | broad (sd 400 bp) | shift | 0.0109 | 0.0109 | 0.0109 |
| mixed 2v2 | broad (sd 400 bp) | shape_any | 0.0054 | 0.0054 | 0.0041 |
| mixed 2v2 | low depth | shape_kernel | 0.0123 | 0.0123 | 0.0129 |
| mixed 2v2 | low depth | shift | 0.0174 | 0.0174 | 0.0169 |
| mixed 2v2 | low depth | shape_any | 0.0152 | 0.0152 | 0.0146 |
| mixed 3v3 | narrow (sd 60 bp) | shape_kernel | 0.0206 | 0.0181 | 0.0164 |
| mixed 3v3 | narrow (sd 60 bp) | shift | 0.0010 | 0.0060 | 0.0054 |
| mixed 3v3 | narrow (sd 60 bp) | shape_any | 0.0109 | 0.0124 | 0.0114 |
| mixed 3v3 | broad (sd 400 bp) | shape_kernel | 0.0000 | 0.0000 | 0.0001 |
| mixed 3v3 | broad (sd 400 bp) | shift | 0.0115 | 0.0117 | 0.0111 |
| mixed 3v3 | broad (sd 400 bp) | shape_any | 0.0053 | 0.0054 | 0.0059 |
| mixed 3v3 | low depth | shape_kernel | 0.0120 | 0.0186 | 0.0161 |
| mixed 3v3 | low depth | shift | 0.0156 | 0.0115 | 0.0144 |
| mixed 3v3 | low depth | shape_any | 0.0136 | 0.0152 | 0.0156 |
| mixed_large 2v2 | narrow (sd 60 bp) | shape_kernel | 0.0233 | 0.0145 | 0.0138 |
| mixed_large 2v2 | narrow (sd 60 bp) | shift | 0.0030 | 0.0077 | 0.0075 |
| mixed_large 2v2 | narrow (sd 60 bp) | shape_any | 0.0136 | 0.0115 | 0.0102 |
| mixed_large 2v2 | broad (sd 400 bp) | shape_kernel | 0.0000 | 0.0089 | 0.0097 |
| mixed_large 2v2 | broad (sd 400 bp) | shift | 0.0114 | 0.0109 | 0.0113 |
| mixed_large 2v2 | broad (sd 400 bp) | shape_any | 0.0053 | 0.0093 | 0.0099 |
| mixed_large 2v2 | low depth | shape_kernel | 0.0125 | 0.0115 | 0.0133 |
| mixed_large 2v2 | low depth | shift | 0.0174 | 0.0132 | 0.0126 |
| mixed_large 2v2 | low depth | shape_any | 0.0160 | 0.0126 | 0.0128 |
| mixed_large 3v3 | narrow (sd 60 bp) | shape_kernel | 0.0207 | 0.0125 | 0.0132 |
| mixed_large 3v3 | narrow (sd 60 bp) | shift | 0.0016 | 0.0097 | 0.0098 |
| mixed_large 3v3 | narrow (sd 60 bp) | shape_any | 0.0115 | 0.0120 | 0.0115 |
| mixed_large 3v3 | broad (sd 400 bp) | shape_kernel | 0.0000 | 0.0096 | 0.0101 |
| mixed_large 3v3 | broad (sd 400 bp) | shift | 0.0127 | 0.0110 | 0.0108 |
| mixed_large 3v3 | broad (sd 400 bp) | shape_any | 0.0056 | 0.0103 | 0.0105 |
| mixed_large 3v3 | low depth | shape_kernel | 0.0116 | 0.0117 | 0.0115 |
| mixed_large 3v3 | low depth | shift | 0.0162 | 0.0117 | 0.0108 |
| mixed_large 3v3 | low depth | shape_any | 0.0146 | 0.0120 | 0.0114 |
<!-- /gen:perclass -->

### 3.4 Region-count boundary

Mixed design (narrow / broad / low depth in equal parts), share of unchanged regions with p < 0.01 over all
classes and tests (nominal 0.01). Sources: `round3/mixed_n_v3`, `round3/mixed_n_v24` (development seeds) and
`round4/confirm_v3_rerun` (seeds 9 100 000+).

| design | regions (total) | p01 range, v2.4 null (dev) | p01 range, v3 null (dev) | p01 range, v3 null (new seeds) |
|---|---|---|---|---|
| 2v2 | 600 | 0.0000–0.0223 | 0.0000–0.0223 | 0.0000–0.0258 |
| 2v2 | 1000 | 0.0000–0.0217 | 0.0000–0.0217 | 0.0000–0.0235 |
| 2v2 | 1500 | 0.0000–0.0234 | 0.0000–0.0234 | 0.0000–0.0230 |
| 2v2 | 2000 | 0.0000–0.0227 | 0.0000–0.0227 | 0.0000–0.0215 |
| 2v2 | 3000 | 0.0000–0.0220 | 0.0075–0.0144 | 0.0080–0.0150 |
| 2v2 | 6000 | 0.0000–0.0234 | 0.0095–0.0134 | 0.0092–0.0132 |
| 3v3 | 600 | 0.0000–0.0192 | 0.0001–0.0159 | 0.0001–0.0181 |
| 3v3 | 1000 | 0.0000–0.0203 | 0.0072–0.0140 | 0.0075–0.0125 |
| 3v3 | 1500 | 0.0000–0.0201 | 0.0088–0.0123 | 0.0095–0.0127 |
| 3v3 | 2000 | 0.0000–0.0195 | 0.0080–0.0132 | 0.0086–0.0127 |
| 3v3 | 3000 | 0.0000–0.0190 | 0.0099–0.0121 | 0.0091–0.0122 |
| 3v3 | 6000 | 0.0000–0.0194 | 0.0094–0.0124 | 0.0092–0.0125 |

The boundary depends on the number of testable regions and on this class mix; it is not a general threshold.

### 3.5 Per-class FDR and stratified BH

<!-- gen:fdrstrata -->
Class FDP (mean over runs) and whole-list FDP [95 % CI], `mixed_large`, by FDR strata:

| design | FDR strata | test | narrow | broad | low depth | whole list |
|---|---|---|---|---|---|---|
| 2v2 | none (global BH) | shape_kernel | 0.025 | 0.222 | 0.090 | 0.045 [0.042, 0.049] / 0.32 |
| 2v2 | none (global BH) | shift | 0.010 | 0.330 | 0.165 | 0.037 [0.032, 0.042] / 0.14 |
| 2v2 | none (global BH) | shape_any | 0.017 | 0.256 | 0.103 | 0.040 [0.036, 0.043] / 0.29 |
| 2v2 | spread | shape_kernel | 0.048 | 0.008 | 0.108 | 0.059 [0.055, 0.062] / 0.34 |
| 2v2 | spread | shift | 0.018 | 0.000 | 0.143 | 0.032 [0.029, 0.036] / 0.16 |
| 2v2 | spread | shape_any | 0.036 | 0.000 | 0.110 | 0.048 [0.044, 0.051] / 0.32 |
| 2v2 | depth | shape_kernel | 0.034 | 0.260 | 0.083 | 0.052 [0.048, 0.056] / 0.31 |
| 2v2 | depth | shift | 0.015 | 0.364 | 0.131 | 0.034 [0.030, 0.038] / 0.14 |
| 2v2 | depth | shape_any | 0.024 | 0.291 | 0.104 | 0.041 [0.038, 0.045] / 0.27 |
| 3v3 | none (global BH) | shape_kernel | 0.024 | 0.211 | 0.060 | 0.045 [0.042, 0.047] / 0.44 |
| 3v3 | none (global BH) | shift | 0.015 | 0.343 | 0.104 | 0.042 [0.038, 0.045] / 0.19 |
| 3v3 | none (global BH) | shape_any | 0.021 | 0.234 | 0.062 | 0.043 [0.040, 0.046] / 0.41 |
| 3v3 | spread | shape_kernel | 0.039 | 0.014 | 0.059 | 0.044 [0.041, 0.046] / 0.43 |
| 3v3 | spread | shift | 0.023 | 0.000 | 0.085 | 0.033 [0.029, 0.036] / 0.19 |
| 3v3 | spread | shape_any | 0.034 | 0.000 | 0.067 | 0.041 [0.038, 0.044] / 0.40 |
| 3v3 | depth | shape_kernel | 0.030 | 0.208 | 0.056 | 0.046 [0.043, 0.048] / 0.43 |
| 3v3 | depth | shift | 0.022 | 0.330 | 0.077 | 0.040 [0.037, 0.044] / 0.18 |
| 3v3 | depth | shape_any | 0.028 | 0.220 | 0.056 | 0.044 [0.042, 0.047] / 0.40 |
<!-- /gen:fdrstrata -->

BH controls the FDR of the whole list. Within-stratum BH (`inference.fdr_strata`) targets the FDR within
strata that match the subgroup of interest; the combined list is then not controlled at the same level (two
strata with one null each: 1 − 0.95² ≈ 0.0975 at 0.05).

## 4. Real data

### 4.1 In-silico insertion on real replicates

Within-condition splits of real replicates (no true difference) with known changes inserted into 20 % of
regions (`benchmarks/profileseg_insertion.py`; output `benchmarks/results/round3/insertion`, run
2026-10-26T21:36:14, regional.py `15a257546ec0aa00`; the shape tests have not changed since).

| data | test | real null: share called / FWER | amplitude-only: share called | shift 20 / 50 / 100 bp: power (FDP) | broaden ×1.3 / ×1.5: power (FDP) |
|---|---|---|---|---|---|
| PRO-seq 3 vs 3 | shape_any | 0.0000 / 0.00 | 0.0000 | 0.94 (0.015) / 1.00 (0.014) / 0.99 (0.015) | 0.23 (0.008) / 0.52 (0.014) |
| PRO-seq 3 vs 3 | shape_kernel | 0.0000 / 0.00 | 0.0000 | 0.36 (0.010) / 0.91 (0.013) / 0.99 (0.014) | 0.12 (0.005) / 0.40 (0.007) |
| PRO-seq 3 vs 3 | shift | 0.0000 / 0.00 | 0.0000 | 0.96 (0.017) / 1.00 (0.016) / 0.35 (0.015) | 0.21 (0.012) / 0.44 (0.015) |
| yeast MNase 2 vs 2 | shape_any | 0.0080 / 0.11 | 0.0050 | 1.00 (0.099) / 0.93 (0.099) / 0.97 (0.101) | 0.98 (0.106) / 0.96 (0.111) |
| yeast MNase 2 vs 2 | shape_kernel | 0.0167 / 0.17 | 0.0160 | 0.10 (0.125) / 0.77 (0.110) / 0.98 (0.111) | 0.98 (0.119) / 0.96 (0.123) |
| yeast MNase 2 vs 2 | shift | 0.0000 / 0.00 | 0.0000 | 1.00 (0.051) / 0.88 (0.051) / 0.34 (0.054) | 0.20 (0.046) / 0.21 (0.045) |

Real-null splits reuse the same replicates; they are not independent validations.

### 4.2 The yeast excess and replicate structure

All excess yeast calls come from the "replicates 1+2 vs 3+4" splits of four genotypes. A label-free measure
(across-pair / within-pair distance of per-sample meta-profiles) tracks them:

| genotype | split | pairing ratio | shape_kernel null calls | shape_any null calls | shift null calls | FDP kernel / any / shift |
|---|---|---|---|---|---|---|
| WT | 1+2 vs rest | 0.96 | 0.00 % | 0.00 % | 0.00 % | 0.081 / 0.059 / 0.050 |
| WT | 1+3 vs rest | 0.93 | 0.00 % | 0.00 % | 0.00 % | 0.026 / 0.046 / 0.055 |
| WT | 1+4 vs rest | 1.13 | 0.00 % | 0.00 % | 0.00 % | 0.047 / 0.033 / 0.025 |
| dKO | 1+2 vs rest | 1.44 | 8.68 % | 5.54 % | 0.00 % | 0.521 / 0.397 / 0.118 |
| dKO | 1+3 vs rest | 0.86 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.015 / 0.040 |
| dKO | 1+4 vs rest | 0.83 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.012 / 0.010 |
| rtf1-8-9A | 1+2 vs rest | 0.96 | 0.00 % | 0.00 % | 0.00 % | 0.033 / 0.060 / 0.063 |
| rtf1-8-9A | 1+3 vs rest | 1.01 | 0.00 % | 0.00 % | 0.00 % | 0.053 / 0.039 / 0.031 |
| rtf1-8-9A | 1+4 vs rest | 1.03 | 0.00 % | 0.00 % | 0.00 % | 0.031 / 0.033 / 0.027 |
| rtf1-11A | 1+2 vs rest | 1.56 | 0.00 % | 0.00 % | 0.00 % | 0.224 / 0.209 / 0.158 |
| rtf1-11A | 1+3 vs rest | 0.82 | 0.00 % | 0.00 % | 0.00 % | 0.008 / 0.019 / 0.030 |
| rtf1-11A | 1+4 vs rest | 0.82 | 0.00 % | 0.00 % | 0.00 % | 0.004 / 0.008 / 0.011 |
| rtf1 | 1+2 vs rest | 1.66 | 3.14 % | 0.00 % | 0.00 % | 0.488 / 0.363 / 0.100 |
| rtf1 | 1+3 vs rest | 0.81 | 0.00 % | 0.00 % | 0.00 % | 0.003 / 0.018 / 0.022 |
| rtf1 | 1+4 vs rest | 0.79 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.006 / 0.014 |
| isw1rtf1 | 1+2 vs rest | 2.28 | 18.28 % | 8.80 % | 0.00 % | 0.602 / 0.515 / 0.097 |
| isw1rtf1 | 1+3 vs rest | 0.72 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.014 / 0.018 |
| isw1rtf1 | 1+4 vs rest | 0.72 | 0.00 % | 0.00 % | 0.00 % | 0.000 / 0.014 / 0.018 |

Splits with pairing ratio > 1.2: 4 of 18; their mean `shape_kernel` FDP 0.459, null calls in 3 of them. Other splits: mean FDP 0.020, null calls in 0 of 14. Spearman (ratio, FDP) = 0.91.

Held-out regions (ratio on one half, tests on the other; predictions fixed beforehand; `round4/holdout`):
both predictions met (null-call share 8.93 % vs 0.00 %;
Spearman 0.51). Interpretation: the excess is associated with a sample-level
difference shared across loci; its origin (batch, culture or other) is unknown — SRA metadata carry no batch
information. The ratio is an aggregate between-group statistic on the same samples, so this is not independent
evidence that the null is valid on these data, and no split is excluded from any result.

## 5. Translation: detection versus size

### 5.1 Power and estimate accuracy (`round4/shift`; seeds 9 400 000+; current decision rule)

| shape | design | d (bp) | power: shift / shape_kernel / shape_any | consistent share of called | direction (consistent) | estimate ÷ true (consistent) | earlier template fit ÷ true (sign ok) |
|---|---|---|---|---|---|---|---|
| gauss60 | 3v3 | 50 | 0.85 / 0.87 / 0.90 | 0.99 | 1.000 | 1.03 | 1.01 (0.90) |
| gauss60 | 3v3 | 100 | 1.00 / 1.00 / 1.00 | 1.00 | 1.000 | 1.00 | 0.83 (0.69) |
| gauss60 | 3v3 | 150 | 0.84 / 1.00 / 1.00 | 1.00 | 1.000 | 0.99 | 0.11 (0.53) |
| gauss60 | 3v3 | 200 | 0.40 / 1.00 / 1.00 | 1.00 | 1.000 | 0.99 | -0.04 (0.50) |
| gauss150 | 3v3 | 50 | 0.07 / 0.17 / 0.15 | 0.97 | 1.000 | 1.38 | 1.37 (0.98) |
| gauss150 | 3v3 | 100 | 0.80 / 0.92 / 0.92 | 0.99 | 1.000 | 1.02 | 1.02 (0.99) |
| gauss150 | 3v3 | 150 | 0.97 / 1.00 / 1.00 | 1.00 | 1.000 | 1.00 | 1.00 (0.97) |
| gauss150 | 3v3 | 200 | 1.00 / 1.00 / 1.00 | 1.00 | 1.000 | 0.99 | 0.97 (0.95) |
| gauss400 | 3v3 | 50 | 0.00 / 0.00 / 0.00 | – | – | – | – (–) |
| gauss400 | 3v3 | 100 | 0.01 / 0.10 / 0.07 | 0.95 | 1.000 | 1.16 | 1.10 (1.00) |
| gauss400 | 3v3 | 150 | 0.17 / 0.60 / 0.55 | 0.97 | 1.000 | 1.04 | 0.99 (1.00) |
| gauss400 | 3v3 | 200 | 0.43 / 0.87 / 0.84 | 0.98 | 1.000 | 1.01 | 0.96 (1.00) |
| bimodal | 3v3 | 50 | 0.94 / 0.80 / 0.94 | 0.98 | 1.000 | 1.01 | 0.84 (0.88) |
| bimodal | 3v3 | 100 | 1.00 / 1.00 / 1.00 | 1.00 | 1.000 | 1.00 | 0.86 (0.80) |
| bimodal | 3v3 | 150 | 0.78 / 1.00 / 1.00 | 1.00 | 1.000 | 0.99 | 0.75 (0.72) |
| bimodal | 3v3 | 200 | 0.07 / 1.00 / 1.00 | 0.06 | 1.000 | 0.96 | 0.27 (0.61) |
| real_GSE112560 | 3v3 | 50 | 0.18 / 0.06 / 0.15 | 0.79 | 1.000 | 1.13 | 1.02 (0.89) |
| real_GSE112560 | 3v3 | 100 | 0.79 / 0.88 / 0.91 | 0.87 | 1.000 | 0.94 | 0.88 (0.87) |
| real_GSE112560 | 3v3 | 150 | 0.82 / 0.99 / 0.98 | 0.93 | 1.000 | 0.94 | 0.87 (0.83) |
| real_GSE112560 | 3v3 | 200 | 0.64 / 0.99 / 0.99 | 0.93 | 1.000 | 0.95 | 0.88 (0.82) |

- Direction of consistent calls: correct in every cell evaluated.
- Estimate (median over consistent calls, d ≥ 50 bp): 0.94–1.08 of the truth in cells with `shift` power ≥ 0.3;
  1.13–1.38 in cells with power ≤ 0.26 (selection among the few regions that pass).
- The `shift` test loses power once the translation exceeds ~1.5–2 peak widths (narrow peaks at 200 bp,
  bimodal at 200 bp); `shape_kernel` and `shape_any` stay monotone. Use `shape_any` for detection.
- The earlier template fit failed beyond ~2 peak widths (sign ≈ 50 %); it is kept only as `shift_bp_fit`.

### 5.2 Not a mechanism: one-sided gains

A peak of half the apex added on one side (not a translation):

| shape | design | offset (bp) | shift power | consistent share of `shift`-called | shape_any power |
|---|---|---|---|---|---|
| gauss60 | 3v3 | 75 | 0.048 | 0.97 | 0.12 |
| gauss60 | 3v3 | 100 | 0.046 | 0.98 | 0.29 |
| gauss60 | 3v3 | 150 | 0.002 | 1.00 | 0.59 |
| gauss150 | 3v3 | 75 | 0.000 | – | 0.01 |
| gauss150 | 3v3 | 100 | 0.001 | 1.00 | 0.02 |
| gauss150 | 3v3 | 150 | 0.001 | 0.60 | 0.09 |
| bimodal | 3v3 | 75 | 0.008 | 0.68 | 0.02 |
| bimodal | 3v3 | 100 | 0.007 | 0.73 | 0.04 |
| bimodal | 3v3 | 150 | 0.000 | – | 0.08 |

`shift` rarely reaches significance for such changes, but when it does most are flagged consistent. A stricter
R² cut (0.8) was chosen on development seeds and rejected on new seeds (`round5/shift_r2_0.8_rejected`): gains
still passed in 5 cells while real-shape translations kept a direction in only 10–70 %.

## 6. Comparison with edgeR diffSpliceDGE

Same simulated counts; edgeR 4.0.16; 50 runs per cell; seeds 9 800 000+ (`round6/diffsplice`, with the
source diff of the run in `src.patch`). Primary endpoint (fixed before the run): end-to-end power at BH 5 %
(every truly changed input region in the denominator); comparator `dS_support_agg60_F` (support bins from the
background-corrected pooled profile, merged to 60 bp, gene-level F), fixed after round 5.

| scenario | design | shape_any | dS_support_agg60_F | difference [paired 95 % CI] | pre-specified hypothesis |
|---|---|---|---|---|---|
| shift50_sd60 | 2v2 | 0.709 (FDP 0.048) | 0.376 (FDP 0.024) | +0.334 [+0.316, +0.353] | met |
| shift50_sd60 | 3v3 | 0.895 (FDP 0.045) | 0.712 (FDP 0.033) | +0.183 [+0.170, +0.196] | met |
| shift100_sd60 | 2v2 | 0.997 (FDP 0.051) | 0.989 (FDP 0.028) | +0.009 [+0.006, +0.011] |  |
| shift100_sd60 | 3v3 | 1.000 (FDP 0.053) | 0.999 (FDP 0.028) | +0.001 [+0.000, +0.001] |  |
| shift50_sd150 | 2v2 | 0.020 (FDP 0.009) | 0.005 (FDP 0.020) | +0.016 [+0.001, +0.032] |  |
| shift50_sd150 | 3v3 | 0.138 (FDP 0.039) | 0.026 (FDP 0.019) | +0.112 [+0.087, +0.135] |  |
| shift100_sd150 | 2v2 | 0.769 (FDP 0.045) | 0.523 (FDP 0.014) | +0.245 [+0.228, +0.264] | met |
| shift100_sd150 | 3v3 | 0.930 (FDP 0.047) | 0.823 (FDP 0.016) | +0.107 [+0.098, +0.116] | met |
| shift200_sd150 | 2v2 | 0.997 (FDP 0.047) | 0.996 (FDP 0.017) | +0.001 [-0.001, +0.002] |  |
| shift200_sd150 | 3v3 | 1.000 (FDP 0.044) | 1.000 (FDP 0.015) | +0.000 [+0.000, +0.000] |  |
| broaden_x1.5 | 2v2 | 0.011 (FDP 0.005) | 0.005 (FDP 0.000) | +0.006 [-0.004, +0.019] |  |
| broaden_x1.5 | 3v3 | 0.224 (FDP 0.031) | 0.068 (FDP 0.020) | +0.156 [+0.134, +0.178] |  |
| side_gain150 | 2v2 | 0.003 (FDP 0.004) | 0.004 (FDP 0.030) | -0.002 [-0.006, +0.005] |  |
| side_gain150 | 3v3 | 0.062 (FDP 0.033) | 0.030 (FDP 0.014) | +0.031 [+0.013, +0.050] |  |

Power–FDP curves (end-to-end power / observed FDP):

| scenario | design | method | BH 0.01 | BH 0.02 | BH 0.05 | BH 0.1 | BH 0.2 |
|---|---|---|---|---|---|---|---|
| shift50_sd60 | 2v2 | shape_any | 0.00 / 0.000 | 0.46 / 0.013 | 0.71 / 0.048 | 0.80 / 0.094 | 0.89 / 0.197 |
| shift50_sd60 | 2v2 | dS_support_agg60_F | 0.16 / 0.007 | 0.23 / 0.012 | 0.38 / 0.024 | 0.50 / 0.052 | 0.64 / 0.104 |
| shift50_sd60 | 3v3 | shape_any | 0.74 / 0.007 | 0.82 / 0.016 | 0.89 / 0.045 | 0.94 / 0.092 | 0.97 / 0.184 |
| shift50_sd60 | 3v3 | dS_support_agg60_F | 0.52 / 0.006 | 0.60 / 0.015 | 0.71 / 0.033 | 0.79 / 0.062 | 0.86 / 0.117 |
| shift100_sd150 | 2v2 | shape_any | 0.01 / 0.000 | 0.60 / 0.013 | 0.77 / 0.045 | 0.85 / 0.101 | 0.91 / 0.204 |
| shift100_sd150 | 2v2 | dS_support_agg60_F | 0.32 / 0.003 | 0.41 / 0.007 | 0.52 / 0.014 | 0.62 / 0.035 | 0.73 / 0.075 |
| shift100_sd150 | 3v3 | shape_any | 0.82 / 0.007 | 0.88 / 0.016 | 0.93 / 0.047 | 0.96 / 0.092 | 0.98 / 0.183 |
| shift100_sd150 | 3v3 | dS_support_agg60_F | 0.69 / 0.002 | 0.75 / 0.004 | 0.82 / 0.016 | 0.88 / 0.035 | 0.92 / 0.077 |
| broaden_x1.5 | 2v2 | shape_any | 0.00 / 0.000 | 0.00 / 0.000 | 0.01 / 0.005 | 0.11 / 0.048 | 0.25 / 0.159 |
| broaden_x1.5 | 2v2 | dS_support_agg60_F | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.01 / 0.004 | 0.03 / 0.013 |
| broaden_x1.5 | 3v3 | shape_any | 0.00 / 0.000 | 0.05 / 0.003 | 0.22 / 0.031 | 0.34 / 0.074 | 0.48 / 0.168 |
| broaden_x1.5 | 3v3 | dS_support_agg60_F | 0.02 / 0.000 | 0.03 / 0.005 | 0.07 / 0.020 | 0.11 / 0.026 | 0.18 / 0.054 |

Controls without a shape change (share of runs with any call):

| control | design | shape_any | dS_all_F | dS_support_F | dS_support_agg60_F |
|---|---|---|---|---|---|
| null | 2v2 | 0.00 | 0.00 | 0.00 | 0.00 |
| null | 3v3 | 0.00 | 0.02 | 0.02 | 0.00 |
| amplitude_x2 | 2v2 | 0.00 | 1.00 | 0.04 | 0.04 |
| amplitude_x2 | 3v3 | 0.00 | 1.00 | 0.02 | 0.00 |
| site_amplitude | 2v2 | 0.00 | 1.00 | 0.00 | 0.04 |
| site_amplitude | 3v3 | 0.00 | 1.00 | 0.04 | 0.00 |

Reading: at BH 5 % diffSplice is conservative (FDP 0.014–0.033) and ProfileSeg runs at 0.045–0.048, so part of
the BH-5 % difference is calibration. At similar observed FDP: 50 bp on narrow peaks — 0.71 vs 0.50 (2 vs 2),
0.89 vs 0.79 (3 vs 3); 100 bp on medium peaks — 0.77 vs 0.73 and 0.93 vs 0.92; broadening at 3 vs 3 — 0.22–0.34
vs 0.18. All-bin diffSplice calls under amplitude-only changes because the raw count distribution over bins
does change when a peak scales over constant background: a difference in target, not an FDR failure. Limits:
ProfileSeg's simulator, comparator configurations chosen by us, no real-count comparison.

## 7. Limitations and open items

1. Independent real data (not used during development), with real-count comparison to diffSplice.
2. The yeast MNase excess: unresolved; replicate-structure diagnostic proposed for the output.
3. Distinguishing translations from other redistributions (R² does not; a model comparison is a candidate).
4. Selection-adjusted or interval estimates for the bp effect.
5. Per-class calibration with few regions (a null that adjusts for the covariate without coarsening the p-value grid).
6. Per-region absolute amplitude test.
7. Spatial dependence beyond duplicated windows; multimodal and asymmetric real shapes beyond the semi-real cells.
8. All simulation evidence comes from one simulator family developed alongside the method.

## 8. Graphical abstract

`docs/figures/profileseg-graphical-abstract-v1.png` (conceptual; prompt in the `.prompt.txt` next to it).
Known issues against the evidence above: the footer's "validation is strongest for selected MNase-seq and
PRO-seq settings" overstates MNase (yeast 2 vs 2 is above nominal); "replicattion" is misspelt; the "Position
shift" panel should be read as an effect size, not as identification of a translation.

## 9. Reproduction

```
python -m pytest -q tests/
PYTHONPATH=benchmarks python benchmarks/profileseg_calibration_boundary.py --out OUT/dev_v3 --semireal SEMIREAL.npz --reps 100
PYTHONPATH=benchmarks python benchmarks/profileseg_calibration_boundary.py --out OUT/confirm_v3 --semireal SEMIREAL.npz --reps 100 --seed-base 9100000
PYTHONPATH=benchmarks python benchmarks/profileseg_calibration_boundary.py --out OUT/mixed_n_v3 --semireal SEMIREAL.npz --reps 100 --factors mixed_n
python benchmarks/profileseg_insertion.py --out OUT/insertion --yeast-cache Y.npz --proseq-cache P.npz
python benchmarks/profileseg_replicate_structure.py Y.npz OUT/insertion/insertion_runs.json OUT/REPLICATE_STRUCTURE.md
python benchmarks/profileseg_replicate_holdout.py Y.npz OUT/holdout
python benchmarks/profileseg_shift_boundary.py --out OUT/shift --templates REAL_TEMPLATES.npz --reps 30 --seed-base 9400000
python benchmarks/profileseg_diffsplice_compare.py --out OUT/diffsplice --reps 50 --seed-base 9800000
python benchmarks/profileseg_review_tables.py benchmarks/results/calibration docs/PROFILESEG_V3_REVIEW.md docs/PROFILESEG_CALIBRATION.md
```

Inputs not in the repository: `SEMIREAL.npz` (GSE112560 bigWigs via `benchmarks/semireal.py`), the yeast
(GSE280671) and PRO-seq (GSE233930) caches (built by `profileseg_insertion.py --build-yeast / --build-proseq`),
`REAL_TEMPLATES.npz` (GSE112560). Results produced before v3.1 record commit `5a5be13` plus an uncommitted diff
identified by its hash; from v3.1 on, results can be tied to a commit.
