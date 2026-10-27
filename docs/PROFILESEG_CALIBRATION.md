# ProfileSeg — calibration boundaries of the per-region shape tests (v3)

Script: `benchmarks/profileseg_calibration_boundary.py`. Run outputs (JSON with provenance, and the
Markdown generated from it): `benchmarks/results/calibration/<run>/`. Every table between `gen:` markers
in this file is generated from those JSON files by
`python benchmarks/profileseg_review_tables.py benchmarks/results/calibration docs/PROFILESEG_CALIBRATION.md`.

## Scope of the claims

- These are **empirical** results in the simulated and semi-real scenarios listed below. "ok" means
  that no statistically clear exceedance of the nominal level was observed in a cell (the 95 % Monte
  Carlo CI of the mean FDP did not lie above 0.05). It does not establish FDR ≤ 0.05.
- The covariate-stratified null (§4) does not depend on the labels, so the exchangeability of label
  assignments within a region is untouched. Comparability of the null across regions is a separate,
  empirical property; it is improved, not guaranteed.

## 1. Protocol

Fixed before the first run on seeds 5 000 000+: parametric shapes, 600 regions, 20 % changed (half
translated by 100 bp, half broadened ×1.5), 2 vs 2 and 3 vs 3, 100 runs per cell. One factor varied at a
time: peak width (sd 60 / 150 / 400 bp), depth (5 / 15 / 40 / 160 reads per bin at the apex), replicate
noise (low / default / high), share of changed regions (0 / 5 / 20 / 50 %), overlapping windows (every
region duplicated 200 bp away), mixed region classes (narrow, broad, low depth; 200 each, reported per
class), a balanced sample-wide 15 bp artefact, real replicate noise (GSE112560). Outcomes at nominal FDR
0.05: mean FDP with 95 % Monte Carlo CI (FWER with an exact binomial CI under the complete null); power;
p01 = share of unchanged regions with p < 0.01 (nominal 0.01).

Deviations from the original protocol: runs per cell raised from 30 to 100 after a one-run smoke test;
`mixed_large` (1000 regions per class) added after the first run; the v3 null was designed after the
first run. Development runs therefore reuse the seeds the design was developed on. The confirmatory run
(§6) uses new seeds and criteria fixed before it was run.

## 2. Runs

<!-- gen:provenance -->
| run | purpose | seeds | null | FDR | commit | uncommitted src files | regional.py sha256 | src diff sha256 | started → finished |
|---|---|---|---|---|---|---|---|---|---|
| confirm_v3 | confirmatory, v3 default (frozen) | 9100000+ | strata | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T21:05:10 → 2026-10-26T21:17:20 |
| dev_fdr_depth | development, depth FDR strata | 5000000+ | strata | depth | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T21:01:48 → 2026-10-26T21:05:09 |
| dev_fdr_spread | development, spread FDR strata | 5000000+ | strata | spread | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:58:23 → 2026-10-26T21:01:46 |
| dev_v24 | development, v2.4 null | 5000000+ | none | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:46:07 → 2026-10-26T20:58:22 |
| dev_v3 | development, v3 default | 5000000+ | strata | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:33:37 → 2026-10-26T20:46:06 |
<!-- /gen:provenance -->

All runs used the same code (identical regional.py hash and identical hash of the uncommitted `src/`
diff on top of commit `5a5be13`).

## 3. Whole-list FDR

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

No cell showed a statistically clear exceedance of 0.05 with either null, in development or
confirmatory runs. The largest point estimates (0.049–0.052) are at the nominal level; the data are
consistent with FDR close to or below 0.05 in these scenarios, not with a demonstrated bound.

## 4. Per-class tail calibration and the v3 null

With the v2.4 null (one null pooled over all regions), p-values of unchanged regions were not
comparable across region types: narrow peaks ≈ 2× liberal in `shape_kernel`, broad peaks fully
conservative, low-depth regions ≈ 1.6× liberal in `shift`. The whole-list FDR stayed near nominal
because conservative and liberal classes offset each other. Cause: F = D / (W + s0) with one global s0 is
compressed where a region's spread is small relative to s0 and heavy-tailed where it is large.

v3 pools the null within quantile strata of each region's total spread (identical under every label
assignment), only when each stratum keeps ≥ 200 regions and ≥ 2000 null values (a coarser p-value grid
otherwise prevents BH from calling small sets). Rejected alternatives: location–scale standardisation
along the covariate; a 2-D grid with a noise-whiteness covariate.

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

Reading:

- **3000 regions (`mixed_large`)**: all classes and tests within 0.0075–0.0145 (0.75–1.45× nominal)
  with the v3 null, in development and confirmatory runs; v2.4 ranged 0.000–0.023.
- **600 regions (`mixed`)**: at 2 vs 2 the v3 null cannot stratify (identical to v2.4): narrow-peak
  `shape_kernel` 0.022–0.026 (≈ 2.2–2.6×), low-depth `shift` ≈ 1.7×. At 3 vs 3 with three strata,
  low-depth `shape_kernel` is 0.016–0.019 and narrow 0.016–0.018. **Small region sets remain
  miscalibrated per class.**
- Low-depth `shift` at 3000 regions: 0.011–0.013 (1.1–1.3×).

## 5. Per-class FDR and stratified BH

Global BH targets the FDR of the whole list (formal control requires valid input p-values and suitable
dependence; here calibration is assessed empirically). A class with almost no power (broad peaks, power ≤ 0.06) collects
mostly false calls even with calibrated p-values. Stratified BH targets the FDR within each stratum, but
(i) only strata that match the subgroup help, and (ii) the combined list is then **not** FDR-controlled
at the same level: two strata with one true null each give 1 − 0.95² ≈ 0.0975 at level 0.05 (the
implementation reproduces 0.099 in 20 000 draws).

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

Spread strata remove the broad-class false calls but leave low-depth regions at 0.06–0.14 and raise
the whole-list `shape_kernel` FDP at 2 vs 2 to 0.059 [0.055, 0.062]; depth strata changed the class FDPs
little in this design (why the depth quantiles did not isolate the low-depth class was not investigated). Decision:
global BH stays the default; `inference.fdr_strata` is an option for a predefined subgroup, with the
combined-list caveat.

## 6. Confirmatory run

Pre-specified (2026-10-26, before the run): frozen code (hashes in §2), new seed base 9 100 000, same
cells, v3 default. C1: no EXCEEDS in any whole-list cell; C2: complete null ≤ 5 of 100 runs with any
call; C3: in `mixed_large`, p01 within [0.007, 0.014] for `shape_kernel` and `shape_any` in every class,
≤ 0.015 for `shift`. Power and per-class FDP reported, not judged.

<!-- gen:confirm -->
| criterion | result |
|---|---|
| C1 whole list: no EXCEEDS | **pass** |
| C2 complete null: ≤ 5 of 100 runs with any call | **pass** |
| C3 per-class tails (mixed_large): p01 in [0.007, 0.014] (shape_kernel, shape_any), ≤ 0.015 (shift) | **pass** |
<!-- /gen:confirm -->

C3 was defined on `mixed_large` only; the 600-region per-class miscalibration (§4) is reproduced in the
confirmatory run and is a stated limitation, not a pass.

## 7. What can and cannot be claimed

- Can: in the scenarios tested (width, depth, noise, share changed, overlapping windows, balanced
  artefacts, real noise; 2 and 3 replicates), whole-list FDR estimates are close to or below 0.05 with no
  clear exceedance, reproduced on new seeds; per-class tails within the pre-specified tolerance at 3000 regions
  (and from 1000 regions at 3 vs 3; at 2 vs 2 only from ~2000 testable regions, see the round-3 sweep in
  `benchmarks/results/round3/`), i.e. within
  ~1.5× of nominal for all three shape tests.
- Cannot: a proven FDR bound; exact permutation validity across regions; FDR within a subgroup (predefined
  or post hoc) under global BH; per-class calibration with a few hundred regions; anything about
  multimodal or asymmetric real shapes beyond the semi-real cell; per-region absolute amplitude, which is
  miscalibrated under complete nulls in two designs (FWER = FDR 0.20 and 0.275).
- Real data: PRO-seq 3 vs 3 is within nominal, but on yeast MNase 2 vs 2 (all 18 within-genotype splits)
  `shape_kernel` has FDP 0.11–0.13 and calls in 17 % of real-null runs (`benchmarks/results/insertion_v3/`).
  A batch explanation for part of this was proposed after seeing the results and is not established.
