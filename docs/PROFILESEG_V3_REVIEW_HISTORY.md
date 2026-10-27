# ProfileSeg v3 — review history (rounds 1–6)

*Complete record of the six review rounds, including withdrawn claims and superseded analyses. The current, consolidated packet is `docs/PROFILESEG_V3_REVIEW.md` (v3.1). Section numbers below are those used during the rounds.*

*2026-10-26. Branch `profileseg-v3`; baseline tag `profileseg-v2.4` (`13206cc`); checkpoint `5a5be13`
plus uncommitted changes, identified by the hashes recorded in every run (§4). Version `3.0.0.dev0`.
"Seg" currently means region-level screening of differential profiles, not segmentation of curves.*

## 0. Round 6 — response to round 5 review

| # | finding | assessment | change / evidence |
|---|---|---|---|
| 1 | bug and export resolved | agreed | — |
| 2 | power denominators differ between methods (filtered regions drop out) | correct in principle | end-to-end power (denominator: every truly changed input region) is now the primary metric; common-set power and regions retained are reported. In practice every method returned ≥ 599 of 600 regions and the two powers differ by ≤ 0.004, so round-5 conclusions are not driven by this |
| 3 | matched-FDP interpolation is post hoc and joins non-monotone branches (9 curves) | correct (all 9 are low-power curves where many runs make no call) | interpolation dropped; pre-specified primary analysis on new seeds (seeds 9 800 000+): end-to-end power at BH 5 %, shape_any versus the round-5-fixed comparator `dS_support_agg60_F`, paired bootstrap CI of the difference; full power–FDP curves with CIs as secondary |
| 4 | advantage limited to some medium translations; 2 vs 2 broadening absolutely weak | correct, and the curves sharpen it (§0a) | at comparable observed FDP the advantage is robust for 50 bp translations of narrow peaks, small for 100 bp on sd 150 at 2 vs 2 and essentially absent at 3 vs 3; 2 vs 2 broadening power at BH 5 % is 0.01 for ProfileSeg |
| 5 | bp is a consistency-checked matching shift, not a mechanism | agreed | wording in docs and report; tiers renamed `high_r2` / `moderate_r2` |
| 6 | round-5 diffSplice run recorded an earlier regional.py hash | correct | the round-6 run stores the full uncommitted source diff (`src.patch`) next to its results; regional.py `87f0c6c86498c6a5` |

### 0a. Pre-specified comparison with edgeR diffSpliceDGE (edgeR 4.0.16; 50 runs per cell; seeds 9 800 000+)

Primary: end-to-end power at BH 5 %; difference shape_any − dS_support_agg60_F with paired bootstrap CI;
hypothesis (fixed before the run): CI > 0 for shift50_sd60 and shift100_sd150, 2 vs 2 and 3 vs 3.

| scenario | design | shape_any | dS_support_agg60_F | difference [paired bootstrap 95 % CI] | hypothesis |
|---|---|---|---|---|---|
| shift20_sd60 | 2v2 | 0.000±0.000 | 0.002±0.001 | -0.002 [-0.003, -0.001] |  |
| shift20_sd60 | 3v3 | 0.001±0.003 | 0.004±0.002 | -0.002 [-0.005, +0.001] |  |
| shift50_sd60 | 2v2 | 0.709±0.017 | 0.376±0.019 | +0.334 [+0.316, +0.353] | **met** |
| shift50_sd60 | 3v3 | 0.895±0.009 | 0.712±0.012 | +0.183 [+0.170, +0.196] | **met** |
| shift100_sd60 | 2v2 | 0.997±0.001 | 0.989±0.002 | +0.009 [+0.006, +0.011] |  |
| shift100_sd60 | 3v3 | 1.000±0.000 | 0.999±0.001 | +0.001 [+0.000, +0.001] |  |
| shift200_sd60 | 2v2 | 1.000±0.000 | 1.000±0.000 | -0.000 [-0.001, +0.000] |  |
| shift200_sd60 | 3v3 | 1.000±0.000 | 1.000±0.000 | -0.000 [-0.001, +0.000] |  |
| shift20_sd150 | 2v2 | 0.000±0.000 | 0.001±0.001 | -0.001 [-0.001, +0.000] |  |
| shift20_sd150 | 3v3 | 0.000±0.000 | 0.000±0.000 | -0.000 [-0.001, +0.000] |  |
| shift50_sd150 | 2v2 | 0.020±0.016 | 0.005±0.002 | +0.016 [+0.001, +0.032] |  |
| shift50_sd150 | 3v3 | 0.138±0.025 | 0.026±0.006 | +0.112 [+0.087, +0.135] |  |
| shift100_sd150 | 2v2 | 0.769±0.017 | 0.523±0.020 | +0.245 [+0.228, +0.264] | **met** |
| shift100_sd150 | 3v3 | 0.930±0.006 | 0.823±0.011 | +0.107 [+0.098, +0.116] | **met** |
| shift200_sd150 | 2v2 | 0.997±0.002 | 0.996±0.001 | +0.001 [-0.001, +0.002] |  |
| shift200_sd150 | 3v3 | 1.000±0.000 | 1.000±0.000 | +0.000 [+0.000, +0.000] |  |
| broaden_x1.5 | 2v2 | 0.011±0.012 | 0.005±0.003 | +0.006 [-0.004, +0.019] |  |
| broaden_x1.5 | 3v3 | 0.224±0.019 | 0.068±0.012 | +0.156 [+0.134, +0.178] |  |
| side_gain150 | 2v2 | 0.003±0.006 | 0.004±0.002 | -0.002 [-0.006, +0.005] |  |
| side_gain150 | 3v3 | 0.062±0.020 | 0.030±0.007 | +0.031 [+0.013, +0.050] |  |

All four pre-specified hypotheses are met. Secondary: end-to-end power / observed FDP along BH levels:

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

Reading. At BH 5 % diffSplice runs at FDP 0.014–0.033 and ProfileSeg at 0.045–0.048, so part of the BH-5 %
difference is calibration, not modelling. Comparing points of similar observed FDP along the curves:
50 bp on sd 60 — ProfileSeg 0.71 (FDP 0.048) vs diffSplice 0.50 (0.052) at 2 vs 2, 0.89 (0.045) vs 0.79 (0.062)
at 3 vs 3: a robust advantage. 100 bp on sd 150 — 0.77 (0.045) vs 0.73 (0.075) at 2 vs 2 and 0.93 (0.047) vs
0.92 (0.077) at 3 vs 3: small to none; the BH-5 % gap there is mostly diffSplice's conservativeness.
Broadening at 3 vs 3: 0.22 (0.031) / 0.34 (0.074) vs 0.18 (0.054): modest; at 2 vs 2 neither method detects
broadening. These curve comparisons use the simulated truth to read FDP and are descriptive. Limits remain:
ProfileSeg's own simulator, comparator configurations chosen by us, no real counts.

## 0‴. Round 5 — response to round 4 review

| # | finding | assessment | change / evidence |
|---|---|---|---|
| 1 | broadening data in the diffSplice comparison rescaled with the inverse library ratio | **correct; our bug** (counts_b × (library_a / library_b)²) | fixed; broadening rerun (§0a). Only the broadening scenario concatenated two simulations |
| 2 | `translation_reliable` cannot mean "a translation happened" | correct | renamed `translation_consistent` (fit-consistency flag); docs and report state that one-sided gains can pass. A stricter R² cut (0.8, chosen on development seeds 7 700 000+) was evaluated on new seeds 9 700 000+ and **rejected** (§0b); R² is now reported as a descriptive tier `translation_fit` (high_r2 ≥ 0.8 / moderate_r2; renamed in round 6 so that it is not read as a confirmed mechanism) |
| 3 | new diagnostics not exported | correct | `region_tests.tsv` now carries `shift_bp_xcorr`, `shift_bp_fit`, `translation_r2`, `translation_ambiguity`, `translation_consistent`, `translation_flag` (direction_mismatch / poor_translation_fit / ambiguous_match / search_edge) and `translation_fit`; unit-tested; report text updated |
| 4 | diffSplice comparison not yet fair | correct; and the round-4 numbers **overstated ProfileSeg's advantage** | support for diffSplice from background-corrected pooled counts (was raw counts); six edgeR configurations; 50 runs per cell with CIs; power compared at matched observed FDP (§0a). The amplitude-only calls of all-bin diffSplice are now described as a difference in target |
| 5 | hold-out check is not independent validation | correct | wording narrowed ("reduces dependence on same-region noise"; region hold-out consistency check; does not resolve the yeast excess) |

Code: regional.py `28584ac7dd480341` (round 4) → current; the decision rule for
direction and bp is unchanged from round 4 (R² ≥ 0.5), so the round-4 translation evaluation (§0′a below,
seeds 9 400 000+) applies; this round adds only exported fields, the failure reason and the R² tier.

### 0‴a. ProfileSeg v3 versus edgeR diffSpliceDGE, corrected (round 5; primary analysis superseded by round 6 §0a) (edgeR 4.0.16; 50 runs per cell; seeds 9 500 000+)

Power at matched observed FDP 0.05 (interpolated along each method's BH-level curve 0.01–0.2; – = not reached). **Superseded in round 6**: the interpolation joins non-monotone curve branches in 9 low-power curves, uses the truth to pick the operating point, and its denominators were method-specific; kept for the record only.
"Best diffSplice" is the best of six configurations per cell, chosen after the fact (favours diffSplice):

| scenario | design | shape_any | shape_kernel | best diffSplice (config) | diffSplice support_F |
|---|---|---|---|---|---|
| shift50_sd60 | 2v2 | 0.70 | 0.62 | 0.51 (support_agg60_F) | 0.36 |
| shift50_sd60 | 3v3 | 0.89 | 0.88 | 0.76 (support_agg60_F) | 0.67 |
| shift100_sd60 | 2v2 | 1.00 | 1.00 | 0.99 (support_agg40_F) | 0.98 |
| shift100_sd60 | 3v3 | 1.00 | 1.00 | 1.00 (support_agg40_F) | 1.00 |
| shift50_sd150 | 2v2 | 0.09 | 0.07 | 0.01 (support_agg40_F) | – |
| shift50_sd150 | 3v3 | 0.21 | 0.22 | 0.10 (support_agg60_F) | – |
| shift100_sd150 | 2v2 | 0.78 | 0.77 | 0.68 (support_agg60_F) | 0.46 |
| shift100_sd150 | 3v3 | 0.93 | 0.93 | 0.88 (support_agg60_F) | 0.76 |
| shift200_sd150 | 2v2 | 0.99 | 1.00 | 1.00 (support_agg60_F) | 0.99 |
| shift200_sd150 | 3v3 | 1.00 | 1.00 | 1.00 (support_agg60_F) | 1.00 |
| broaden_x1.5 | 2v2 | 0.09 | 0.12 | 0.04 (support_agg60_F) | – |
| broaden_x1.5 | 3v3 | 0.26 | 0.34 | 0.18 (support_agg60_F) | 0.07 |
| side_gain150 | 2v2 | 0.05 | 0.06 | 0.01 (support_agg60_F) | 0.00 |
| side_gain150 | 3v3 | 0.12 | 0.16 | 0.05 (support_agg60_F) | 0.01 |

At BH 5 % (power ± 95 % CI / FDP ± CI):

| scenario | design | shape_any | dS_support_F | dS_support_agg60_F | dS_all_F |
|---|---|---|---|---|---|
| shift50_sd60 | 2v2 | 0.69±0.02 / 0.044±0.009 | 0.26±0.02 / 0.021±0.007 | 0.38±0.02 / 0.024±0.006 | 0.01±0.00 / 0.017±0.023 |
| shift50_sd60 | 3v3 | 0.89±0.01 / 0.048±0.006 | 0.60±0.02 / 0.028±0.005 | 0.71±0.02 / 0.032±0.005 | 0.07±0.01 / 0.005±0.006 |
| shift100_sd150 | 2v2 | 0.78±0.02 / 0.050±0.007 | 0.31±0.02 / 0.011±0.005 | 0.52±0.02 / 0.013±0.004 | 0.10±0.01 / 0.018±0.012 |
| shift100_sd150 | 3v3 | 0.93±0.01 / 0.046±0.006 | 0.66±0.01 / 0.014±0.004 | 0.82±0.01 / 0.020±0.004 | 0.45±0.02 / 0.020±0.004 |
| broaden_x1.5 | 2v2 | 0.02±0.02 / 0.006±0.007 | 0.00±0.00 / 0.000±0.000 | 0.01±0.00 / 0.000±0.000 | 0.00±0.00 / 0.040±0.055 |
| broaden_x1.5 | 3v3 | 0.24±0.03 / 0.041±0.011 | 0.03±0.01 / 0.007±0.010 | 0.07±0.01 / 0.016±0.011 | 0.04±0.01 / 0.033±0.024 |

Controls without a shape change (share of runs with any call):

| control | design | shape_any | dS_all_F | dS_support_F | dS_support_agg60_F |
|---|---|---|---|---|---|
| null | 2v2 | 0.00 | 0.04 | 0.02 | 0.02 |
| null | 3v3 | 0.00 | 0.00 | 0.00 | 0.00 |
| amplitude_x2 | 2v2 | 0.00 | 1.00 | 0.00 | 0.02 |
| amplitude_x2 | 3v3 | 0.00 | 1.00 | 0.02 | 0.02 |
| site_amplitude | 2v2 | 0.00 | 1.00 | 0.00 | 0.04 |
| site_amplitude | 3v3 | 0.00 | 1.00 | 0.06 | 0.00 |

Reading:
- With a background-corrected support and bin aggregation, diffSplice is much stronger than in round 4
  (50 bp on sd 60 at 3 vs 3: 0.30 → 0.60 support_F, 0.71 with 60 bp bins at BH 5 %). The round-4 gap was
  largely a configuration artefact.
- At matched FDP 0.05 ProfileSeg remains more sensitive for medium translations (e.g. 0.89 vs 0.76 and
  0.70 vs 0.51 for 50 bp on sd 60; 0.93 vs 0.88 and 0.78 vs 0.68 for 100 bp on sd 150) and for broadening
  (0.26, `shape_kernel` 0.34, vs 0.18 at 3 vs 3); large translations are ≈ 1.00 for both.
- diffSplice remains conservative at BH 5 % (FDP 0.01–0.03); ProfileSeg's largest FDP is 0.055 ± 0.007.
- Limits: ProfileSeg's own simulator; edgeR configurations chosen by us (not exhaustive); no real counts yet.
  Supportable claim: more sensitive in these configurations and simulations, not generally superior.

### 0‴b. Stricter translation cut evaluated and rejected (seeds 9 700 000+)

R² ≥ 0.8 (development seeds: one-sided gains 0.15 kept, translations 0.75) on new seeds: direction among
consistent calls 1.000 everywhere (pass); estimate bias outside [0.85, 1.15] in 7 cells (was 4; the cut does not
reduce selection bias); one-sided gains above 0.2 in 5 cells (gauss60 35 bp 3v3: 0.37; gauss60 50 bp 3v3: 0.52; gauss60 75 bp 3v3: 0.38; gauss60 100 bp 3v3: 0.23; gauss150 200 bp 3v3: 0.23) — **fail**; and translations
kept with a direction:

| shape | 3v3 kept share at d = 50 / 100 / 150 / 200 bp |
|---|---|
| gauss60 | 0.65 / 0.97 / 1.00 / 1.00 |
| gauss150 | 0.48 / 0.63 / 0.89 / 0.97 |
| gauss400 | – / 0.07 / 0.30 / 0.45 |
| bimodal | 0.40 / 0.97 / 1.00 / 0.07 |
| real_GSE112560 | 0.10 / 0.31 / 0.57 / 0.70 |

On real shapes most translations would lose their bp value while gains still pass in places. R² does not
separate the two reliably, so it is reported (`translation_r2`, `translation_fit`) rather than used as a decision
rule. Distinguishing translations from other redistributions remains open.

## 0′. Round 4 — response to round 3 review

| # | finding | assessment | change / evidence |
|---|---|---|---|
| 1 | translation estimate fails at large shifts; direction near random; label and bp value can disagree | **correct for the bp value; the Higher/Lower label was right**. The label comes from the projection sign, which was correct in 100 % of called translated regions in every cell; the per-sample template fit (a local search from zero) put the sign right in only ~50–60 % at 150–200 bp on narrow peaks | **method change**: `estimate` is now the cross-correlation shift of the condition means (all lags, no capture range); Higher/Lower only when the translation is reliable (same sign as the projection, `translation_r2` ≥ 0.5, `translation_ambiguity` < 0.9), else `Changed` with NaN bp; the template fit is kept as `shift_bp_fit`. Evaluated against criteria fixed beforehand on new seeds (§0a): direction **pass**; bias **fails** at the detection limit (selection); one-sided gains **fail** (flagged reliable when `shift` is significant) |
| 2 | pairing ratio uses the same curves and the evaluated split; cannot declare those splits invalid nulls | **correct; claim withdrawn**. The ratio is an aggregate between-group statistic on the same data, so its correlation with the calls is expected | wording: "excess calls are associated with replicate structure; cause undetermined". A held-out-region check (ratio on half A, tests on half B, predictions written beforehand) met both predictions (§0b); this reduces the dependence on noise in the same regions (regions may still be correlated across halves) — the samples are shared |
| 3 | `shape_any` is not a universally reliable final list | correct | stated: on yeast it inherits the excess of the structured splits; a combination cannot repair a miscalibrated component |
| 4 | region boundary should be read as testable regions, specific mix | correct | the output now reports testable regions, number of null strata and null values per stratum (note in the per-region results); the boundary was reproduced on new seeds (seeds 9 100 000+: 2 vs 2 per-class p01 0.000–0.022 up to 2000 regions, 0.008–0.015 at 3000; 3 vs 3 0.0075–0.0125 from 1000) and is described as specific to this mix |
| 5 | edgeR diffSplice comparison missing for v3 | correct | same simulated counts, same regions (§0c) |

Code changed this round (regional.py `28584ac7dd480341`). The confirmatory calibration run was repeated
with this code on its original seeds: all 50 original cells identical (`benchmarks/results/round4/confirm_v3_rerun`).

### 0′a. Translation estimate (criteria fixed before the run; seeds 9 400 000+)

Criteria: (1) direction accuracy among reliable called translations ≥ 0.95; (2) median reliable estimate ÷ true
within [0.85, 1.15] for d ≥ 50 bp; (3) one-sided gain: ≤ 0.2 of `shift`-called regions flagged reliable.
Source: `benchmarks/results/round4/shift` (all shapes, both designs, d = 10–200 bp).

| shape | design | d (bp) | `shift` power | reliable share of called | direction of reported estimate (reliable) | reported estimate ÷ true (reliable) | template fit ÷ true (sign ok) |
|---|---|---|---|---|---|---|---|
| gauss60 | 3v3 | 50 | 0.85 | 0.99 | 1.000 | 1.03 | 1.01 (0.90) |
| gauss60 | 3v3 | 100 | 1.00 | 1.00 | 1.000 | 1.00 | 0.83 (0.69) |
| gauss60 | 3v3 | 150 | 0.84 | 1.00 | 1.000 | 0.99 | 0.11 (0.53) |
| gauss60 | 3v3 | 200 | 0.40 | 1.00 | 1.000 | 0.99 | -0.04 (0.50) |
| bimodal | 3v3 | 50 | 0.94 | 0.98 | 1.000 | 1.01 | 0.84 (0.88) |
| bimodal | 3v3 | 100 | 1.00 | 1.00 | 1.000 | 1.00 | 0.86 (0.80) |
| bimodal | 3v3 | 150 | 0.78 | 1.00 | 1.000 | 0.99 | 0.75 (0.72) |
| bimodal | 3v3 | 200 | 0.07 | 0.06 | 1.000 | 0.96 | 0.27 (0.61) |
| real_GSE112560 | 3v3 | 50 | 0.18 | 0.79 | 1.000 | 1.13 | 1.02 (0.89) |
| real_GSE112560 | 3v3 | 100 | 0.79 | 0.87 | 1.000 | 0.94 | 0.88 (0.87) |
| real_GSE112560 | 3v3 | 150 | 0.82 | 0.93 | 1.000 | 0.94 | 0.87 (0.83) |
| real_GSE112560 | 3v3 | 200 | 0.64 | 0.93 | 1.000 | 0.95 | 0.88 (0.82) |

- (1) **pass** in every cell (1.000).
- (2) **fails in 4 cells**, all where `shift` power ≤ 0.07 (estimate ÷ true 1.16–1.38): selection bias among the
  few regions that pass the threshold. Elsewhere 0.94–1.03 for unimodal and real shapes.
- (3) **fails**: one-sided gains rarely reach significance in `shift` (≤ 5 % of changed regions), but when they
  do, 60–100 % are flagged reliable with an estimate of 0.3–0.8 × the offset. `translation_r2` does not separate
  a one-sided gain from a translation. Not re-tuned on these seeds.

| shape | design | offset (bp) | `shift` power (share of changed regions called) | reliable share of called | `shape_any` power |
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

### 0′b. Held-out-region check of the replicate-structure association

| split | pairing ratio (A) | null calls on B: shape_kernel / shape_any / shift |
|---|---|---|
| isw1rtf1:isw1rtf1_1+isw1rtf1_2 | 1.90 | 23.06 % / 12.12 % / 0.00 % |
| rtf1:rtf1_1+rtf1_2 | 1.51 | 3.44 % / 3.18 % / 0.00 % |
| dKO:dKO_1+dKO_2 | 1.41 | 9.23 % / 6.15 % / 0.00 % |
| rtf1-11A:rtf1-11A_1+rtf1-11A_2 | 1.35 | 0.00 % / 0.00 % / 0.00 % |
| WT:WT_1+WT_4 | 1.09 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-8-9A:rtf1-8-9A_1+rtf1-8-9A_3 | 1.00 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-8-9A:rtf1-8-9A_1+rtf1-8-9A_2 | 1.00 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-8-9A:rtf1-8-9A_1+rtf1-8-9A_4 | 1.00 | 0.00 % / 0.00 % / 0.00 % |
| WT:WT_1+WT_3 | 0.96 | 0.00 % / 0.00 % / 0.00 % |
| WT:WT_1+WT_2 | 0.96 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-11A:rtf1-11A_1+rtf1-11A_4 | 0.89 | 0.00 % / 0.00 % / 0.00 % |
| dKO:dKO_1+dKO_3 | 0.87 | 0.00 % / 0.00 % / 0.00 % |
| rtf1-11A:rtf1-11A_1+rtf1-11A_3 | 0.85 | 0.00 % / 0.00 % / 0.00 % |
| rtf1:rtf1_1+rtf1_3 | 0.85 | 0.00 % / 0.00 % / 0.00 % |
| dKO:dKO_1+dKO_4 | 0.84 | 0.00 % / 0.00 % / 0.00 % |
| rtf1:rtf1_1+rtf1_4 | 0.81 | 0.00 % / 0.00 % / 0.00 % |
| isw1rtf1:isw1rtf1_1+isw1rtf1_3 | 0.77 | 0.00 % / 0.00 % / 0.00 % |
| isw1rtf1:isw1rtf1_1+isw1rtf1_4 | 0.76 | 0.00 % / 0.00 % / 0.00 % |

Pre-specified predictions: (1) higher null-call share on B for ratio(A) > 1.2: **met** (8.93 % vs 0.00 %; splits with calls 3 of 4 vs 0 of 14); (2) Spearman > 0.5: **met** (0.51).


Reading: the association holds when the ratio and the calls come from disjoint regions; this is a region hold-out consistency check, not validation on independent samples, and it does not resolve the excess error rate on yeast (rtf1-11A 1+2 has ratio
1.35 and no calls). The same replicates enter both halves, so this shows a sample-level difference shared
across loci, of unknown origin; it is not evidence that the pooled null is valid on these data, and the splits
are not excluded from any result. Spearman 0.51 only just meets the pre-specified 0.5.

### 0′c. ProfileSeg v3 versus edgeR QL diffSpliceDGE (round 4; superseded by §0a — broadening data wrong, support on raw counts)

| scenario | design | shape_any | shape_kernel | shift | diffSplice_all | diffSplice_support | ProfileSeg translation |
|---|---|---|---|---|---|---|---|
| shift20_sd60 | 2v2 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.050 | 0.00 · nan · nan |
| shift20_sd60 | 3v3 | 0.00 / 0.000 | 0.00 / 0.000 | 0.01 / 0.017 | 0.00 / 0.000 | 0.00 / 0.000 | 0.19 · 1.91 · 1.00 |
| shift50_sd60 | 2v2 | 0.69 / 0.045 | 0.59 / 0.040 | 0.63 / 0.049 | 0.01 / 0.017 | 0.06 / 0.027 | 0.99 · 1.09 · 1.00 |
| shift50_sd60 | 3v3 | 0.89 / 0.048 | 0.86 / 0.036 | 0.84 / 0.047 | 0.08 / 0.006 | 0.30 / 0.021 | 0.99 · 1.03 · 1.00 |
| shift100_sd60 | 2v2 | 1.00 / 0.043 | 0.99 / 0.043 | 0.97 / 0.040 | 0.72 / 0.021 | 0.89 / 0.028 | 1.00 · 0.99 · 1.00 |
| shift100_sd60 | 3v3 | 1.00 / 0.057 | 1.00 / 0.052 | 0.99 / 0.056 | 0.96 / 0.021 | 0.99 / 0.028 | 1.00 · 0.99 · 1.00 |
| shift200_sd60 | 2v2 | 1.00 / 0.052 | 1.00 / 0.049 | 0.32 / 0.063 | 1.00 / 0.020 | 1.00 / 0.031 | 1.00 · 0.98 · 1.00 |
| shift200_sd60 | 3v3 | 1.00 / 0.051 | 1.00 / 0.054 | 0.40 / 0.035 | 1.00 / 0.019 | 1.00 / 0.030 | 1.00 · 0.99 · 1.00 |
| shift20_sd150 | 2v2 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.050 | 0.00 / 0.000 | 0.00 · nan · nan |
| shift20_sd150 | 3v3 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.04 · 2.66 · 1.00 |
| shift50_sd150 | 2v2 | 0.03 / 0.019 | 0.06 / 0.034 | 0.01 / 0.014 | 0.00 / 0.000 | 0.00 / 0.000 | 0.10 · 1.53 · 1.00 |
| shift50_sd150 | 3v3 | 0.14 / 0.034 | 0.16 / 0.024 | 0.06 / 0.033 | 0.00 / 0.000 | 0.01 / 0.000 | 0.59 · 1.41 · 1.00 |
| shift100_sd150 | 2v2 | 0.78 / 0.046 | 0.77 / 0.054 | 0.54 / 0.041 | 0.09 / 0.022 | 0.30 / 0.022 | 0.99 · 1.05 · 1.00 |
| shift100_sd150 | 3v3 | 0.93 / 0.046 | 0.92 / 0.040 | 0.82 / 0.043 | 0.44 / 0.019 | 0.68 / 0.021 | 0.99 · 1.02 · 1.00 |
| shift200_sd150 | 2v2 | 0.99 / 0.055 | 0.99 / 0.045 | 0.94 / 0.056 | 0.98 / 0.021 | 0.99 / 0.019 | 1.00 · 1.00 · 1.00 |
| shift200_sd150 | 3v3 | 1.00 / 0.050 | 1.00 / 0.053 | 1.00 / 0.047 | 1.00 / 0.019 | 1.00 / 0.019 | 1.00 · 1.00 · 1.00 |
| broaden_x1.5 | 2v2 | 0.03 / 0.013 | 0.11 / 0.051 | 0.00 / 0.000 | 0.00 / 0.050 | 0.01 / 0.000 | – |
| broaden_x1.5 | 3v3 | 0.25 / 0.034 | 0.35 / 0.043 | 0.00 / 0.000 | 0.02 / 0.027 | 0.07 / 0.013 | – |
| side_gain150 | 2v2 | 0.01 / 0.013 | 0.05 / 0.037 | 0.00 / 0.014 | 0.00 / 0.000 | 0.00 / 0.050 | – |
| side_gain150 | 3v3 | 0.08 / 0.014 | 0.13 / 0.037 | 0.00 / 0.000 | 0.00 / 0.000 | 0.01 / 0.000 | – |
| null | 2v2 | 0.00 | 0.00 | 0.00 | 0.10 | 0.05 | – |
| null | 3v3 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | – |
| amplitude_x2 | 2v2 | 0.00 | 0.00 | 0.00 | 1.00 | 0.05 | – |
| amplitude_x2 | 3v3 | 0.00 | 0.00 | 0.00 | 1.00 | 0.05 | – |
| site_amplitude | 2v2 | 0.00 | 0.00 | 0.00 | 1.00 | 0.00 | – |
| site_amplitude | 3v3 | 0.00 | 0.00 | 0.00 | 1.00 | 0.10 | – |

Reading:
- Translations of 50–100 bp and broadening: ProfileSeg detects more (e.g. 50 bp on sd 60: `shape_any`
  0.69 / 0.89 vs diffSplice-support 0.06 / 0.30; ×1.5 broadening at 3 vs 3: 0.25 vs 0.07). Large translations:
  both ≈ 1.00.
- diffSplice over all core bins calls in every run under amplitude-only changes: with constant background, a
  higher peak does change the raw count distribution over bins, so this is a difference in target, not an FDR
  failure (amended in round 5). Restricted to support bins it does not call. ProfileSeg made no calls in any control.
- diffSplice is conservative here (FDP 0.02–0.03) where ProfileSeg runs at 0.04–0.057; part of the power gap may
  be calibration rather than modelling.
- diffSplice gives no bp estimate; ProfileSeg's reliable translations are within 0.98–1.09 of the truth except
  near the detection limit.
- Fairness limits: data from ProfileSeg's own simulator (Gaussian or template shapes, negative-binomial counts);
  one diffSplice configuration (20 bp bins, no filtering or bin aggregation, gene-level F); no real-data counts
  comparison yet. This is a same-data comparison, not evidence of general superiority.

## 0″. Response to round 2 review (kept for the record; claims amended in round 4)

| # | finding | assessment | evidence (this round) |
|---|---|---|---|
| 1 | "≥ ~1000 regions" overstated: `mixed_large` had 3000 regions; at 2 vs 2 strata need ~2000 testable regions | **correct** | region-count sweep (§6a): at 2 vs 2 per-class calibration is unchanged from v2.4 up to 2000 regions (fewer than 2000 testable, no strata) and improves at 3000 and 6000; at 3 vs 3 it is within 0.7–1.4× from 1000 regions, partial at 600. Claims now state these numbers |
| 2 | real data fail on yeast | **result correct; interpretation partly rebutted**. All excess calls come from the "replicates 1+2 vs 3+4" splits of four genotypes; the 14 other splits have no null calls and FDP 0.000–0.081. A label-free measure that uses no ProfileSeg test (distance between per-sample meta-profiles across vs within the 1,2 / 3,4 pairs) is 1.44–2.28 in exactly those four genotypes and 0.96 in the two without calls; Spearman (ratio, FDP) = 0.91 over 18 splits (§7a). ~~In those splits the two groups differ systematically, so they are not valid nulls.~~ (withdrawn in round 4: the association is not independent evidence) What we cannot show: the origin (batch, culture day or other); SRA metadata carry no batch information (one submission, one instrument model). The ratio and its 1.2 cut-off were examined after the results and are reported as such | 
| 3 | `shape_any` had no real-data validation | **correct** | added to the insertion benchmark (§7): PRO-seq no null calls, FDP ≤ 0.015, power ≥ the better component in most cells; yeast FDP 0.099–0.111 with all splits, inheriting `shape_kernel`'s behaviour on the structured splits (0.006–0.060 elsewhere) |
| 4 | confirmatory scope limited; 1.4–1.5× is a tolerance, not exact calibration | correct | wording: "within the pre-specified tolerance", not "calibrated to nominal"; the confirmatory run shares the simulator, stated as a limitation |
| 5 | `shift` loses power for larger shifts | **correct, and the estimate is worse than noted** (§6b): `shift` power falls once d exceeds ~1.5–2× the peak sd (sd 60 bp: 1.00 at 100 bp, 0.39 at 200 bp; bimodal: 0.06 at 200 bp), while `shape_kernel` and `shape_any` stay monotone (0.99–1.00). The translation estimate collapses in the same range (estimate ÷ true 0.12 and −0.08 at 150 / 200 bp for sd 60 bp; 0.43 for bimodal at 200 bp) and is biased upwards among called regions near the detection limit (1.4–1.9×; selection bias). Partial rebuttal: this is a property of the `shift` test and the estimate, not of shape detection, which `shape_any` covers | 

## 1. Response to round 1 review

| # | finding | assessment | change |
|---|---|---|---|
| 1 | stratified BH docstring claimed "FDR within every stratum and hence overall" | **correct; our claim was wrong**. Two strata with one true null each give 1 − 0.95² ≈ 0.0975 at level 0.05; the implementation gives 0.099 in 20 000 draws | docstring, user docs, config comment and plan corrected: within-stratum FDR only; the combined list is not FDR-controlled at the same level. Global BH remains the default |
| 2 | "controlled in every cell" overstated a no-exceedance criterion | correct | all documents now say "no statistically clear exceedance in the scenarios tested; FDR estimates close to or below 0.05"; the verdict legend states that "ok" does not establish FDR ≤ 0.05. A confirmatory run on new seeds with pre-specified criteria was added (§5) |
| 3 | label invariance of the covariate does not make regions exchangeable | correct | code and docs now separate the two: within-region exchangeability of label assignments is untouched; cross-region comparability is an empirical property, improved but not guaranteed. Low-depth `shift` excess is listed as a limitation, not as solved by the depth warning |
| 4 | numbers in the packet did not match the full v3 table | **correct; our error**. Three "v3" rows (width 400 bp, depth 5, high noise; 3 vs 3) had been copied from the v2.4-null table | every table in this packet and in `PROFILESEG_CALIBRATION.md` is now generated from run JSON files by `benchmarks/profileseg_review_tables.py`; each run records commit, uncommitted-diff hash, regional.py hash, options, seeds and times (§4) |
| 5 | ACAT is not an unconditional guarantee | correct | `shape_any` is described as a Cauchy combination that is approximately valid through a tail approximation, with calibration checked by simulation (it is judged separately in the confirmatory run) |
| 6 | absolute amplitude FWER 0.20 / 0.275 called "mildly liberal" | correct: under a complete null FDR = FWER, so this is clear miscalibration; the share of false calls among exact-null regions does not offset it | wording changed; per-region amplitude inference is not claimed until fixed, and the paper scope is the shape tests |

## 2. What the reviewer is asked to check now

Round 3: whether the region-count boundary (§6a), the shift boundary (§6b), the `shape_any` real-data
results (§7) and the replicate-structure argument (§7a) are read correctly, and whether the argument in
§7a is strong enough to treat the yeast 1+2 vs 3+4 splits as invalid nulls.

Round 2 questions (still open):

1. Whether the corrected statements (§1) and the claims in §8 are supported by the generated tables.
2. Whether the confirmatory run (§5) is adequate: criteria fixed beforehand, frozen code, new seeds;
   whether criterion C3 (defined on 3000 regions, 1000 per class) is a fair choice given the known small-set failure.
3. Whether reporting per-class calibration and per-class FDR separately (§6) is the right framing.
4. Whether any remaining statement overreaches.

## 3. Method in one paragraph

ProfileSeg compares signal profiles around predefined anchors between two conditions with biological
replicates. Shape tests use background-corrected, mass-normalised profiles: `shape_kernel` (an MMD-type
kernel distance between condition means, scaled by within-condition spread plus a label-invariant s0),
`shift` (a directional projection on the pooled template's derivative; effect size is a translation-model
estimate in bp) and `shape_any` (per-region Cauchy combination of the two p-values, then one BH). With 2
or 3 replicates per group there are only 3 or 10 label assignments, so per-region p-values come from a
null pooled over regions × non-observed assignments. In v3 that pooling is done within quantile strata
of each region's total spread (identical under every assignment) when ≥ 200 regions and ≥ 2000 null
values per stratum are available. BH runs over all regions.

## 4. Runs and provenance

<!-- gen:provenance -->
| run | purpose | seeds | null | FDR | commit | uncommitted src files | regional.py sha256 | src diff sha256 | started → finished |
|---|---|---|---|---|---|---|---|---|---|
| confirm_v3 | confirmatory, v3 default (frozen) | 9100000+ | strata | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T21:05:10 → 2026-10-26T21:17:20 |
| dev_fdr_depth | development, depth FDR strata | 5000000+ | strata | depth | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T21:01:48 → 2026-10-26T21:05:09 |
| dev_fdr_spread | development, spread FDR strata | 5000000+ | strata | spread | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:58:23 → 2026-10-26T21:01:46 |
| dev_v24 | development, v2.4 null | 5000000+ | none | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:46:07 → 2026-10-26T20:58:22 |
| dev_v3 | development, v3 default | 5000000+ | strata | none | `5a5be1371e` | 4 | `15a257546ec0aa00` | `2995aee9d62eb521` | 2026-10-26T20:33:37 → 2026-10-26T20:46:06 |
<!-- /gen:provenance -->

Outputs: `benchmarks/results/calibration/<run>/calibration_boundary.json` and the Markdown generated in
the same run. Development runs reuse seeds 5 000 000+, on which the v3 null was designed; they are
development evidence only.

## 5. Whole-list FDR and the confirmatory run

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

Confirmatory criteria, fixed before the run (`docs/PROFILESEG_CALIBRATION.md` §6): C1 no EXCEEDS in any
whole-list cell; C2 complete null ≤ 5 of 100 runs with any call; C3 in `mixed_large`, p01 within
[0.007, 0.014] for `shape_kernel` and `shape_any` in every class, ≤ 0.015 for `shift`.

<!-- gen:confirm -->
| criterion | result |
|---|---|
| C1 whole list: no EXCEEDS | **pass** |
| C2 complete null: ≤ 5 of 100 runs with any call | **pass** |
| C3 per-class tails (mixed_large): p01 in [0.007, 0.014] (shape_kernel, shape_any), ≤ 0.015 (shift) | **pass** |
<!-- /gen:confirm -->

Interpretation: in the scenarios tested, whole-list FDR estimates are close to or below 0.05 with no
clear exceedance, reproduced on new seeds. This is empirical support, not a proven bound.

## 6. Per-class calibration and per-class FDR

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

- With 3000 regions (`mixed_large`) the v3 null brings every class and test to 0.75–1.45× nominal (v2.4: 0–2.3×); the region-count boundary is in §6a.
- With 600 regions it cannot help at 2 vs 2 (no strata; identical to v2.4: narrow `shape_kernel`
  2.2–2.6×) and only partly at 3 vs 3 (low-depth and narrow `shape_kernel` 1.6–1.9×). **Small region sets
  remain miscalibrated per class**; this reproduces in the confirmatory run.

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

Calibrated p-values do not give per-class FDR: the broad class (power ≤ 0.06) has class FDP 0.21–0.34
under global BH. Spread strata fix that class and raise the whole-list FDP (0.059 at 2 vs 2 for
`shape_kernel`); depth strata change little here. A user who will read a predefined low-power subgroup on
its own needs strata matching that subgroup and gives up the combined-list guarantee.

### 6a. Region-count boundary (review round 2, point 1)

Mixed design (narrow / broad / low depth, one third each) at 600–6000 regions in total, 100 runs per
cell, development seeds; share of unchanged regions with p < 0.01 (nominal 0.01), v2.4 null → v3 null.
Source: `benchmarks/results/round3/mixed_n_v3`, `mixed_n_v24`.

| design | regions (total) | test | narrow | broad | low depth | whole-list FDP [CI] / power (v3) |
|---|---|---|---|---|---|---|
| 2v2 | 600 | shape_kernel | 0.022 → 0.022 | 0.000 → 0.000 | 0.012 → 0.012 | 0.044 [0.037, 0.052] / 0.33 |
| 2v2 | 600 | shift | 0.002 → 0.002 | 0.012 → 0.012 | 0.017 → 0.017 | 0.034 [0.025, 0.043] / 0.15 |
| 2v2 | 600 | shape_any | 0.012 → 0.012 | 0.005 → 0.005 | 0.016 → 0.016 | 0.035 [0.028, 0.042] / 0.30 |
| 2v2 | 1000 | shape_kernel | 0.022 → 0.022 | 0.000 → 0.000 | 0.013 → 0.013 | 0.048 [0.042, 0.054] / 0.34 |
| 2v2 | 1000 | shift | 0.003 → 0.003 | 0.011 → 0.011 | 0.016 → 0.016 | 0.035 [0.028, 0.042] / 0.16 |
| 2v2 | 1000 | shape_any | 0.013 → 0.013 | 0.005 → 0.005 | 0.016 → 0.016 | 0.042 [0.036, 0.048] / 0.32 |
| 2v2 | 1500 | shape_kernel | 0.023 → 0.023 | 0.000 → 0.000 | 0.013 → 0.013 | 0.048 [0.042, 0.053] / 0.34 |
| 2v2 | 1500 | shift | 0.003 → 0.003 | 0.011 → 0.011 | 0.018 → 0.018 | 0.037 [0.031, 0.043] / 0.17 |
| 2v2 | 1500 | shape_any | 0.013 → 0.013 | 0.005 → 0.005 | 0.016 → 0.016 | 0.041 [0.036, 0.045] / 0.32 |
| 2v2 | 2000 | shape_kernel | 0.023 → 0.023 | 0.000 → 0.000 | 0.013 → 0.013 | 0.046 [0.041, 0.050] / 0.34 |
| 2v2 | 2000 | shift | 0.003 → 0.003 | 0.011 → 0.011 | 0.016 → 0.016 | 0.039 [0.034, 0.044] / 0.16 |
| 2v2 | 2000 | shape_any | 0.013 → 0.013 | 0.005 → 0.005 | 0.016 → 0.016 | 0.041 [0.037, 0.045] / 0.33 |
| 2v2 | 3000 | shape_kernel | 0.022 → 0.014 | 0.000 → 0.008 | 0.013 → 0.012 | 0.041 [0.038, 0.044] / 0.32 |
| 2v2 | 3000 | shift | 0.003 → 0.008 | 0.011 → 0.011 | 0.017 → 0.013 | 0.035 [0.029, 0.040] / 0.13 |
| 2v2 | 3000 | shape_any | 0.013 → 0.011 | 0.005 → 0.009 | 0.015 → 0.012 | 0.037 [0.033, 0.041] / 0.28 |
| 2v2 | 6000 | shape_kernel | 0.023 → 0.013 | 0.000 → 0.010 | 0.013 → 0.012 | 0.043 [0.041, 0.045] / 0.31 |
| 2v2 | 6000 | shift | 0.003 → 0.010 | 0.011 → 0.011 | 0.017 → 0.012 | 0.037 [0.034, 0.041] / 0.12 |
| 2v2 | 6000 | shape_any | 0.013 → 0.012 | 0.005 → 0.010 | 0.016 → 0.012 | 0.040 [0.037, 0.042] / 0.27 |
| 3v3 | 600 | shape_kernel | 0.019 → 0.016 | 0.000 → 0.000 | 0.010 → 0.015 | 0.040 [0.035, 0.045] / 0.42 |
| 3v3 | 600 | shift | 0.002 → 0.006 | 0.010 → 0.010 | 0.015 → 0.013 | 0.034 [0.026, 0.041] / 0.22 |
| 3v3 | 600 | shape_any | 0.010 → 0.011 | 0.004 → 0.004 | 0.013 → 0.014 | 0.040 [0.035, 0.045] / 0.41 |
| 3v3 | 1000 | shape_kernel | 0.020 → 0.014 | 0.000 → 0.007 | 0.012 → 0.012 | 0.041 [0.037, 0.045] / 0.43 |
| 3v3 | 1000 | shift | 0.001 → 0.008 | 0.012 → 0.012 | 0.016 → 0.012 | 0.044 [0.037, 0.051] / 0.20 |
| 3v3 | 1000 | shape_any | 0.011 → 0.012 | 0.006 → 0.010 | 0.015 → 0.012 | 0.043 [0.038, 0.048] / 0.41 |
| 3v3 | 1500 | shape_kernel | 0.020 → 0.012 | 0.000 → 0.010 | 0.011 → 0.011 | 0.048 [0.045, 0.052] / 0.45 |
| 3v3 | 1500 | shift | 0.002 → 0.009 | 0.012 → 0.010 | 0.015 → 0.012 | 0.040 [0.035, 0.045] / 0.19 |
| 3v3 | 1500 | shape_any | 0.011 → 0.012 | 0.006 → 0.011 | 0.015 → 0.012 | 0.044 [0.040, 0.047] / 0.42 |
| 3v3 | 2000 | shape_kernel | 0.019 → 0.012 | 0.000 → 0.008 | 0.012 → 0.013 | 0.045 [0.041, 0.048] / 0.44 |
| 3v3 | 2000 | shift | 0.001 → 0.009 | 0.011 → 0.011 | 0.017 → 0.012 | 0.039 [0.035, 0.043] / 0.19 |
| 3v3 | 2000 | shape_any | 0.011 → 0.011 | 0.006 → 0.009 | 0.016 → 0.013 | 0.042 [0.039, 0.045] / 0.41 |
| 3v3 | 3000 | shape_kernel | 0.019 → 0.012 | 0.000 → 0.010 | 0.012 → 0.012 | 0.044 [0.042, 0.047] / 0.44 |
| 3v3 | 3000 | shift | 0.002 → 0.010 | 0.012 → 0.011 | 0.015 → 0.011 | 0.040 [0.037, 0.044] / 0.19 |
| 3v3 | 3000 | shape_any | 0.011 → 0.011 | 0.005 → 0.011 | 0.014 → 0.012 | 0.042 [0.040, 0.045] / 0.41 |
| 3v3 | 6000 | shape_kernel | 0.019 → 0.012 | 0.000 → 0.010 | 0.011 → 0.011 | 0.045 [0.043, 0.046] / 0.45 |
| 3v3 | 6000 | shift | 0.002 → 0.009 | 0.012 → 0.011 | 0.016 → 0.011 | 0.042 [0.039, 0.045] / 0.19 |
| 3v3 | 6000 | shape_any | 0.011 → 0.012 | 0.006 → 0.010 | 0.014 → 0.011 | 0.044 [0.042, 0.046] / 0.42 |

At 2 vs 2 nothing changes up to 2000 regions (fewer than 2000 testable: no strata); at 3 vs 3 the
improvement starts at 1000. Whole-list FDP stays at or below 0.048 (CI upper bound ≤ 0.054) everywhere.

### 6b. Shift size × peak shape (review round 2, point 5)

20 % of 600 regions translated by d bp; 30 runs per cell; seeds 9 300 000+. Source:
`benchmarks/results/round3/shift` (all d from 10 to 200 bp there; FDP ≤ 0.064 in every cell).

| shape | design | d (bp) | power: shift / shape_kernel / shape_any | estimate ÷ true: all translated · called |
|---|---|---|---|---|
| gauss60 | 2v2 | 20 | 0.00 / 0.00 / 0.00 | 0.98 · 2.08 |
| gauss60 | 2v2 | 50 | 0.66 / 0.66 / 0.73 | 0.97 · 1.05 |
| gauss60 | 2v2 | 100 | 0.97 / 1.00 / 1.00 | 0.88 · 0.88 |
| gauss60 | 2v2 | 150 | 0.78 / 1.00 / 1.00 | 0.18 · 0.32 |
| gauss60 | 2v2 | 200 | 0.32 / 1.00 / 1.00 | -0.10 · 0.14 |
| gauss60 | 3v3 | 20 | 0.01 / 0.01 / 0.01 | 1.00 · 1.85 |
| gauss60 | 3v3 | 50 | 0.85 / 0.87 / 0.89 | 0.97 · 1.01 |
| gauss60 | 3v3 | 100 | 1.00 / 1.00 / 1.00 | 0.82 · 0.82 |
| gauss60 | 3v3 | 150 | 0.85 / 1.00 / 1.00 | 0.12 · 0.15 |
| gauss60 | 3v3 | 200 | 0.39 / 1.00 / 1.00 | -0.08 · 0.06 |
| gauss150 | 2v2 | 20 | 0.00 / 0.00 / 0.00 | 0.94 · nan |
| gauss150 | 2v2 | 50 | 0.02 / 0.05 / 0.02 | 0.96 · 1.31 |
| gauss150 | 2v2 | 100 | 0.59 / 0.76 / 0.78 | 0.97 · 1.04 |
| gauss150 | 2v2 | 150 | 0.89 / 0.97 / 0.98 | 0.98 · 0.99 |
| gauss150 | 2v2 | 200 | 0.95 / 0.99 / 1.00 | 0.97 · 0.97 |
| gauss150 | 3v3 | 20 | 0.00 / 0.00 / 0.00 | 0.95 · nan |
| gauss150 | 3v3 | 50 | 0.05 / 0.17 / 0.14 | 0.96 · 1.38 |
| gauss150 | 3v3 | 100 | 0.81 / 0.92 / 0.93 | 0.99 · 1.01 |
| gauss150 | 3v3 | 150 | 0.97 / 1.00 / 1.00 | 0.99 · 0.99 |
| gauss150 | 3v3 | 200 | 0.99 / 1.00 / 1.00 | 0.97 · 0.97 |
| gauss400 | 2v2 | 20 | 0.00 / 0.00 / 0.00 | 0.56 · nan |
| gauss400 | 2v2 | 50 | 0.00 / 0.00 / 0.00 | 0.62 · nan |
| gauss400 | 2v2 | 100 | 0.00 / 0.02 / 0.00 | 0.66 · nan |
| gauss400 | 2v2 | 150 | 0.04 / 0.30 / 0.23 | 0.72 · 1.01 |
| gauss400 | 2v2 | 200 | 0.16 / 0.68 / 0.63 | 0.78 · 0.94 |
| gauss400 | 3v3 | 20 | 0.00 / 0.00 / 0.00 | 0.71 · nan |
| gauss400 | 3v3 | 50 | 0.00 / 0.00 / 0.00 | 0.75 · nan |
| gauss400 | 3v3 | 100 | 0.02 / 0.11 / 0.09 | 0.78 · 1.01 |
| gauss400 | 3v3 | 150 | 0.18 / 0.57 / 0.53 | 0.82 · 1.00 |
| gauss400 | 3v3 | 200 | 0.40 / 0.87 / 0.83 | 0.86 · 0.94 |
| bimodal | 2v2 | 20 | 0.01 / 0.00 / 0.00 | 0.78 · 1.42 |
| bimodal | 2v2 | 50 | 0.78 / 0.38 / 0.77 | 0.81 · 0.85 |
| bimodal | 2v2 | 100 | 0.99 / 1.00 / 1.00 | 0.86 · 0.86 |
| bimodal | 2v2 | 150 | 0.64 / 1.00 / 1.00 | 0.82 · 0.80 |
| bimodal | 2v2 | 200 | 0.02 / 1.00 / 1.00 | 0.36 · 0.39 |
| bimodal | 3v3 | 20 | 0.01 / 0.00 / 0.00 | 0.77 · 0.90 |
| bimodal | 3v3 | 50 | 0.93 / 0.79 / 0.94 | 0.83 · 0.84 |
| bimodal | 3v3 | 100 | 1.00 / 1.00 / 1.00 | 0.87 · 0.87 |
| bimodal | 3v3 | 150 | 0.78 / 1.00 / 1.00 | 0.77 · 0.78 |
| bimodal | 3v3 | 200 | 0.06 / 1.00 / 1.00 | 0.43 · 0.40 |
| real_GSE112560 | 2v2 | 20 | 0.00 / 0.00 / 0.00 | 0.69 · nan |
| real_GSE112560 | 2v2 | 50 | 0.04 / 0.01 / 0.02 | 0.77 · 1.06 |
| real_GSE112560 | 2v2 | 100 | 0.59 / 0.61 / 0.72 | 0.83 · 0.89 |
| real_GSE112560 | 2v2 | 150 | 0.68 / 0.91 / 0.93 | 0.87 · 0.90 |
| real_GSE112560 | 2v2 | 200 | 0.48 / 0.97 / 0.96 | 0.86 · 0.91 |
| real_GSE112560 | 3v3 | 20 | 0.00 / 0.00 / 0.00 | 0.78 · nan |
| real_GSE112560 | 3v3 | 50 | 0.16 / 0.10 / 0.16 | 0.81 · 1.00 |
| real_GSE112560 | 3v3 | 100 | 0.80 / 0.89 / 0.91 | 0.85 · 0.87 |
| real_GSE112560 | 3v3 | 150 | 0.83 / 0.99 / 0.99 | 0.88 · 0.90 |
| real_GSE112560 | 3v3 | 200 | 0.62 / 0.99 / 0.99 | 0.81 · 0.88 |

## 7. Real data with the v3 null

In-silico insertion on real replicates (`benchmarks/profileseg_insertion.py`; 20 % of regions changed,
BH 5 %), rerun with the code identified in §4; output and provenance in
`benchmarks/results/insertion_v3/`.

| data | test | null: share called / FWER | amplitude-only: share called | shift 20 bp | shift 50 bp | shift 100 bp | broaden ×1.3 | broaden ×1.5 |
|---|---|---|---|---|---|---|---|---|
| yeast 2 vs 2 | shape_kernel | 0.0167 / 0.17 | 0.0160 | 0.10 / 0.125 | 0.77 / 0.110 | 0.98 / 0.111 | 0.98 / 0.119 | 0.96 / 0.123 |
| yeast 2 vs 2 | shift | 0.0000 / 0.00 | 0.0000 | 1.00 / 0.051 | 0.88 / 0.051 | 0.34 / 0.054 | 0.20 / 0.046 | 0.21 / 0.045 |
| yeast 2 vs 2 | shape_any | 0.0080 / 0.11 | 0.0050 | 1.00 / 0.099 | 0.93 / 0.099 | 0.97 / 0.101 | 0.98 / 0.106 | 0.96 / 0.111 |
| PRO-seq 3 vs 3 | shape_kernel | 0.0000 / 0.00 | 0.0000 | 0.36 / 0.010 | 0.91 / 0.013 | 0.99 / 0.014 | 0.12 / 0.005 | 0.40 / 0.007 |
| PRO-seq 3 vs 3 | shift | 0.0000 / 0.00 | 0.0000 | 0.96 / 0.017 | 1.00 / 0.016 | 0.35 / 0.015 | 0.21 / 0.012 | 0.44 / 0.015 |
| PRO-seq 3 vs 3 | shape_any | 0.0000 / 0.00 | 0.0000 | 0.94 / 0.015 | 1.00 / 0.014 | 0.99 / 0.015 | 0.23 / 0.008 | 0.52 / 0.014 |

Source: `benchmarks/results/round3/insertion` (provenance in `insertion.json`; per-run records in
`insertion_runs.json`). PRO-seq 3 vs 3: no calls under real nulls, FDP ≤ 0.017 for all three tests.
Yeast 2 vs 2 over all 18 splits: `shape_kernel` FDP 0.110–0.125, `shape_any` 0.099–0.111, `shift`
0.045–0.054. Real-null splits reuse the same replicates and are not independent validations.

### 7a. Where the yeast excess comes from

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


Reading (amended in round 4): the excess calls are associated with replicate structure; the cause is undetermined.
The ratio is an aggregate between-group statistic on the same data, so the association is not independent
evidence about the null (see §0b for a held-out-region check); but the ratio and the 1.2 cut-off were chosen
after seeing the results, the origin of the pair structure is unknown, and a user with such structure
and no batch annotation would get these calls. A label-free replicate-structure diagnostic in the
output is a proposed addition (§9).

## 8. Claims

| claim | status |
|---|---|
| whole-list FDR near or below 0.05 across width, depth, noise, share changed, overlap, balanced artefact, real noise (2 and 3 replicates) | supported empirically, reproduced on new seeds; not a proven bound |
| per-class tail calibration within the pre-specified tolerance (≤ ~1.5×) | supported at 3000 regions (2 vs 2 and 3 vs 3) and from 1000 regions at 3 vs 3; at 2 vs 2 not below ~2000 testable regions (§6a) |
| per-class calibration with a few hundred regions | **not supported** (2.2–2.6× liberal for narrow peaks at 2 vs 2) |
| FDR within a subgroup under global BH | **not supported**; stratified BH targets it within matching strata and loses the combined-list level |
| exact permutation validity across regions | not claimed |
| `shape_any` calibration | supported in simulation and on PRO-seq; on yeast it inherits the excess of the structured splits; ACAT validity is approximate |
| per-region absolute amplitude | **miscalibrated** under complete nulls in two designs (FWER = FDR 0.20, 0.275); not claimed |
| translation estimate (bp; cross-correlation since round 4; a consistency-checked matching shift between the two mean profiles, not a mechanism) | for consistent calls on translated regions: direction correct in every evaluated cell; median estimate ÷ true 0.94–1.03 away from the detection limit; biased upwards (1.16–1.5×) near it (selection); the earlier template fit failed beyond ~2 peak widths |
| `translation_consistent` identifies translations | **not supported**: a fit-consistency flag; one-sided gains that reach significance pass it; a stricter R² cut was evaluated and rejected (§0b) |
| more sensitive than edgeR diffSplice | pre-specified primary analysis met for 50 bp (sd 60) and 100 bp (sd 150) translations at BH 5 %; at comparable observed FDP the advantage is robust for 50 bp on narrow peaks, small or absent for 100 bp on sd 150, modest for broadening at 3 vs 3, and neither detects broadening at 2 vs 2; ProfileSeg's simulator and our comparator choice; no real counts |
| `shift` power increases with shift size | **not supported**: falls beyond ~1.5–2× peak sd; use `shape_any` for detection |
| real-data calibration | PRO-seq 3 vs 3 within nominal; yeast 2 vs 2 above nominal over all splits (`shape_kernel` 0.11–0.13, `shape_any` 0.10–0.11), associated with replicate structure, also with held-out regions; cause undetermined (§7a, §0′b) |

## 9. Open items

00. From rounds 4–5: separating translations from one-sided gains needs a different discriminator than R²;
    real-count diffSplice comparison (BAM path).
    From round 4: separate one-sided gains from translations (e.g. compare a translation model with a gain
    model per region); a selection-adjusted or interval estimate for called regions; diffSplice on real counts
    (BAM path) and with tuned configurations; independent real data.

0. Proposed from this round: a label-free replicate-structure diagnostic (within-group pairing ratio)
   reported with the per-region tests; a reliability flag for the translation estimate (d beyond the
   template's support, multimodal ambiguity) and a note on selection bias among called regions.

1. Uncertainty for the shift estimate (interval, fit quality, a "not a translation" flag).
2. Same-target comparisons: csaw / DiffBind (amplitude), MMDiff2 (shape, version and configuration
   stated), IWTomics with aligned inputs (localisation), compared per question.
3. Validation on real data not used during development.
4. Per-class calibration for small region sets (e.g. a smoothed covariate-adjusted null that does not
   coarsen the p-value grid), the residual low-depth `shift` excess, and the absolute amplitude test.
5. Spatial dependence beyond duplicated windows; real multimodal or asymmetric shapes.

## 10. Reproduction

```
python -m pytest -q tests/
PYTHONPATH=benchmarks python benchmarks/profileseg_calibration_boundary.py --out OUT/dev_v3 \
    --semireal SEMIREAL.npz --reps 100                       # also: --null none; --fdr-strata spread|depth
PYTHONPATH=benchmarks python benchmarks/profileseg_calibration_boundary.py --out OUT/confirm_v3 \
    --semireal SEMIREAL.npz --reps 100 --seed-base 9100000
python benchmarks/profileseg_review_tables.py OUT docs/PROFILESEG_V3_REVIEW.md docs/PROFILESEG_CALIBRATION.md
# round 3
PYTHONPATH=benchmarks python benchmarks/profileseg_calibration_boundary.py --out OUT/mixed_n_v3 \
    --semireal SEMIREAL.npz --reps 100 --factors mixed_n          # and --null none
python benchmarks/profileseg_shift_boundary.py --out OUT/shift --templates REAL_TEMPLATES.npz --reps 30
python benchmarks/profileseg_insertion.py --out OUT/insertion --yeast-cache Y.npz --proseq-cache P.npz
python benchmarks/profileseg_replicate_structure.py Y.npz OUT/insertion/insertion_runs.json OUT/REPLICATE_STRUCTURE.md
```

`SEMIREAL.npz` is built from GSE112560 bigWigs by `python benchmarks/semireal.py CONFIG OUT.npz`. Raw
data are not in the repository.
