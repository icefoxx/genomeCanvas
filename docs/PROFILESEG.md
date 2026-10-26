# ProfileSeg — amplitude and shape comparison of aligned profiles

ProfileSeg compares aligned signal profiles (ChIP-seq, CUT&RUN, ATAC-seq) between two conditions. It
keeps **amplitude** (how much signal) separate from **shape** (where and how spread the signal is),
and makes **replicate-level** inferences. It is a research prototype.

* **Global (v1.1):** over a fixed region set, does the typical profile change in amplitude or shape?
* **Per region (v2):** which regions change? Empirical-Bayes moderated tests with BH FDR.
* Not yet: segmentation / local "where along the profile" tests, BAM count models (see
  [PROFILESEG_PLAN.md](PROFILESEG_PLAN.md)). Related work: [PROFILESEG_LITERATURE.md](PROFILESEG_LITERATURE.md).

## Quick start

```bash
pip install -e ".[studio]"                       # GenomeCanvas + pyBigWig/scipy
gcanvas-profileseg example --out ps_example --benchmark-reps 30
```

This single command:
1. writes synthetic bigWigs (3 WT + 3 KO replicates, +100 bp shift), a BED file and a manifest;
2. runs the full analysis, reading the bigWigs through the same extraction code as real data;
3. runs a 30-repetition benchmark with known truth.

Look at `ps_example/results/reports/` and `ps_example/results/figures/`.

### Your own data

```bash
gcanvas-profileseg init --out config.yaml          # annotated template
# edit config.yaml and metadata/sample_manifest.tsv, then:
gcanvas-profileseg run --config config.yaml
```

**`metadata/sample_manifest.tsv`** (tab-separated; leave unknown fields empty and they are
flagged, not guessed):

```
sample_id  condition  biological_replicate  assay  genome_build  signal_file  signal_type  normalization  region_set  batch  pairing_id  bigwig_params  fragment_length
WT_1       WT         1                     ChIP   mm10          bw/WT_1.bw   coverage     CPM            TSS         1      p1          binSize=10;extendReads  210
KO_1       KO         1                     ChIP   mm10          bw/KO_1.bw   coverage     CPM            TSS         1      p1          binSize=10;extendReads  205
...
```

**Regions** must be one fixed, condition-independent BED for all samples. Strand is needed for
oriented anchors (TSS/TES); the `summit` anchor reads narrowPeak column 10. Instead of bigWigs you
can pass `matrix:` pointing at a deepTools `computeMatrix reference-point` file or a ProfileSeg
`.npz`.

### Power guidance

```bash
gcanvas-profileseg power --out ps_power                                   # parametric peaks
gcanvas-profileseg power --templates-config config.yaml --regions 300 1000 # your own shapes
```

Simulates replicates per condition × number of regions × effect size (shift, amplitude fold,
broadening, site-specific fold) and reports power together with the false-call rate of the same design
(`reports/06_power.md`, `results/metrics/power_grid.tsv`). Use a pilot of your own data for realistic
signal-to-noise.

### Semi-synthetic benchmark on your own profile shapes

```bash
gcanvas-profileseg benchmark --templates-config config.yaml --reps 40 --out ps_benchmark
```

The strongest reference-condition profiles (ranked by background-corrected signal, lightly
smoothed) are used as templates. Known changes are injected and replicate noise is simulated on top.
Use this to check power and calibration for *your* signal-to-noise level before trusting a result.

## What is measured

| family | metric | meaning |
|---|---|---|
| amplitude | `log2_auc`, `log2_centre_auc` | log2 fold change of mean AUC (signal × bp) over the window and over ±`center_half` |
| shape | `centroid` | shift of the signal's centre of mass (bp; oriented anchors only) |
| shape | `sd_width`, `iqr_width` | change in spread (bp) |
| shape | `asymmetry` | change in (right − left)/total mass around the anchor (oriented only) |
| shape | `group_w1` | Wasserstein-1 distance (bp) between condition-mean shapes; a distance, so check its permutation p |
| naive | `naive_sd_width`, `naive_centroid` | the same metrics **without** background correction, reported only to show the confound |

**Shape pipeline** (each step exists because the benchmark showed it is needed):
1. **Background:** the *mean* of the outer `background.width` bp on each side, subtracted per profile.
   The median is biased low on count noise.
2. **Core window:** ±`core_flank`.
3. **Inclusion is condition-independent:** a region is kept if its pooled (all-sample) corrected mass
   is positive and above the `min_mass_quantile`.
4. **Bounded scaling:** each region is divided by *pooled mass + median pooled mass*, a factor that
   does not depend on condition.
5. **Aggregation per replicate**, with **signed** moments. No per-region truncation, because clipping
   negative noise makes taller peaks look narrower.

**Inference** is replicate-level:
- each replicate is summarised across regions, then the conditions are contrasted;
- primary interval: **Welch t** on the replicate summaries;
- secondary interval: **hierarchical bootstrap** (replicates within condition, then regions, keeping
  regions paired across conditions);
- **exact permutation** of replicate labels, reporting the **smallest attainable p**. For 3 vs 3 this
  is 0.1, so permutation alone can never reach 0.05; for 1 vs 1 every result is labelled exploratory.

### Roles and multiplicity (v1.1)

* **Primary** metrics (pre-specified: `primary_amplitude`, default `log2_auc`; `primary_shape`, default
  `auto` = centroid for oriented anchors, otherwise iqr_width) are judged on their own Welch interval.
* **Secondary** metrics are Holm-adjusted together; a secondary interval that excludes 0 but fails Holm is
  reported as "no clear change after Holm".
* **Naive** metrics are comparators only and never "called".
* `log2_enrichment` (centre signal / local background) is insensitive to global rescaling; compare it with
  `log2_auc` when the background flag fires.

**Diagnostics** (02_QC.md and flags): background level ratio between conditions; region-level
heterogeneity of amplitude changes (re-weighting caveat); region-level shifts vs a within-condition
replicate-split noise floor (bidirectional shifts look like broadening in the aggregate); replicates
on the wrong side of a called primary effect; replicate concordance matrix.

**Region weighting** (default `shrunk_sample`): each sample's region profile is divided by its own mass +
that sample's median mass. Uniform amplitude changes cancel, and region-specific gains re-weight regions
less than with pooled scaling (iqr false calls 0.85 → 0.15 when the broadest 20 % of regions gain ×2.5).

## Per-region tests (v2)

For every region and metric, replicates are the observations; residual variances are shrunk towards a
prior estimated across regions (Smyth 2004; identical to `limma::eBayes`, checked in the unit tests);
BH FDR per metric (`inference.region_fdr`, default 0.05).

| metric | meaning |
|---|---|
| `log2_auc` | window AUC (absolute; inherits upstream normalisation) |
| `log2_signal` | background-corrected core signal (absolute) |
| `log2_signal_relative` | `log2_signal` minus a per-sample offset estimated on presumed-unchanged regions: **site-specific** changes. Assumes most regions do not change |
| `log2_enrichment` | centre signal / own background |
| `shape_kernel` | **main per-region shape test (v2.1)**: kernel distance between condition-mean profiles, null pooled over regions and replicate-label assignments; Changed / Not_significant, no direction. Works for 2 vs 2 |
| `shift` | per-region shift (bp) from a template fit (oriented anchors): the effect size for position changes. Attenuated by ~10–28 % in validation; do not use for calling at 2 vs 2 on broad profiles |
| `log2_scale` | per-region width scale from the template fit; exploratory |
| `centroid` | per-region centre of mass (oriented anchors); low power on broad profiles |
| `mad_width` | per-region spread (mean absolute deviation); **low power, exploratory** |

Shape tests use every region with positive pooled background-corrected signal (condition-blind).
**Replicate diagnostics:** each replicate's sample-wide shift and width offset relative to its own group
(flags a deviant replicate or pair), and, with ≥ 3 replicates per group, the share of shape calls that
survive leaving each replicate out. Real-data checks (yeast MNase nulls and remodeler mutants; human PRO-seq
NELF degradation) are in `docs/PROFILESEG_COMPARISON.md` §6; application scenarios and a tool comparison in
`docs/PROFILESEG_APPLICATIONS.md`.
Validation of the v2.1 shape tests (frozen code, fresh seeds, real replicate noise, 3 vs 3 and
2 vs 2) is in `docs/PROFILESEG_COMPARISON.md` §4: `shape_kernel` passed every pre-specified
criterion; `shift` failed two 2 vs 2 cells. Two amplitude design decisions, both forced by
benchmark failures:

* **Relative amplitude uses background-corrected signal, not window AUC.** Window AUC includes
  background, so when IP efficiency drops, weak (background-dominated) regions fall less and look
  relatively "up". On GSE112560 TCF3 this produced a strength-dependent pattern (Spearman −0.33) that
  disappeared with background-corrected signal (+0.07). The report prints this correlation.
* **Offsets use the half-sample mode, then iterative exclusion of changed regions.** A plain per-sample
  median is pulled by one-sided changes (FDR 0.30 when 20 % of regions gain); the robust offset matches an
  oracle that knows the changed regions (FDR 0.025).

**Robust global shape:** for each replicate, the median of per-region sd / centroid; Welch across
replicates. Unaffected by re-weighting; when it disagrees with the aggregate width/centroid, a flag says so.

Outputs: `reports/05_regional.md`, `results/metrics/region_tests.tsv` (region × metric: estimate,
moderated t, df, p, FDR, call), `region_tests_summary.tsv`, `robust_global_shape.tsv`,
`figures/G_region_volcano.*`.

## Outputs

```
results/
  reports/00_discovery.md            inputs, contrast, anchor/orientation, metadata flags
  reports/01_methods_and_assumptions.md  formulas, parameters, results (observed / inferred / hypothesised)
  reports/02_QC.md                   inclusion counts, missingness, zero/negative bins, heterogeneity
  reports/03_benchmark.md            (benchmark/example) call rates, bias, CI coverage with denominators
  reports/04_limitations.md
  reports/05_regional.md             per-region tests (v2)
  reports/03b_benchmark_regional.md  (benchmark) per-region FDR and power
  results/metrics/effects.tsv        one row per metric: estimate, Welch CI, bootstrap CI, permutation p and floor
  results/metrics/shape_per_replicate.tsv, amplitude_per_profile.tsv, region_effects.tsv, group_distributions.tsv
  results/profiles.npz               the extracted tensor (reusable as `matrix:` input)
  figures/A_…E_… (.pdf/.svg/.png)    replicate profiles, raw vs normalised, difference curve, effect sizes
  run_log.json                       config, versions, seeds, runtime, flags
```

## Benchmark results (v1 pilot)

40 repetitions per scenario, 3 vs 3 replicates, 300 regions, ±2 kb at 20 bp. Each entry is
**call rate | truth | mean estimate | 95 % CI coverage**. Truth is the same analysis applied to
the noise-free version of each run's regions. Real templates are the 2,000 strongest ZNF384-WT HA
ChIP TSS profiles (GSE112560, mm10), a low signal-to-noise case.

| scenario | metric | truly changes | parametric | real templates |
|---|---|---|---|---|
| null | sd_width | no | 0.03 \| 0 \| −0.31 \| 0.97 | 0.05 \| 0 \| 0.67 \| 0.95 |
| amplitude (peak only) | log2_auc | yes | 1.00 \| 0.224 \| 0.225 \| 1.00 | 1.00 \| 0.319 \| 0.317 \| 0.88 |
| amplitude (peak only) | **sd_width** | **no** | **0.00** \| 0 \| −0.26 \| 1.00 | **0.05** \| 0 \| 0.08 \| 0.95 |
| amplitude (peak only) | naive_sd_width | (confound) | 1.00 \| −71.1 \| −71.1 \| 1.00 | 1.00 \| −32.6 \| −32.4 \| 0.97 |
| global scaling | sd_width | no | 0.05 \| 0 \| 0.98 \| 0.95 | 0.00 \| 0 \| −0.48 \| 1.00 |
| shift +100 bp | centroid | yes | 1.00 \| 100 \| 100 \| 0.95 | 0.95 \| 27.6 \| 28.4 \| 0.97 |
| broaden | iqr_width | yes | 1.00 \| 103 \| 103 \| 0.97 | 0.12 \| −6.6 \| −6.2 \| 0.95 |
| right flank | asymmetry | yes | 1.00 \| 0.281 \| 0.281 \| 0.97 | 1.00 \| 0.252 \| 0.253 \| 0.95 |
| left flank | asymmetry | yes | 1.00 \| −0.281 \| −0.280 \| 0.95 | 1.00 \| −0.267 \| −0.267 \| 0.95 |
| 25 % regions shifted | centroid | yes | 1.00 \| 24.7 \| 24.5 \| 1.00 | 0.10 \| 7.0 \| 6.3 \| 1.00 |
| fragment length (technical) | iqr_width | yes | 0.72 \| 7.8 \| 7.6 \| 1.00 | 0.00 \| −0.3 \| −0.5 \| 1.00 |

**Across all metrics:**
- false-call rate where the truth is ≈ 0: median 0.03 (maximum 0.20 parametric, 0.17 real);
- CI coverage: minimum 0.90 (parametric) and 0.88 (real).

**Read the real-template column carefully.** On weak, broad TSS profiles the *true* effects
themselves are small: a 100 bp peak translation moves the TSS-window centroid by only 28 bp. Power
for effects under about 7 bp is low. Run the template benchmark on your own data to see what is
detectable.

## Limitations (v1)

- **Global target only.** v1 does not say which regions change or where within the profile.
  Per-region effects are descriptive.
- **Few replicates.** With 3 vs 3, Welch intervals rest on few summaries and the permutation floor is
  0.1. With 1 vs 1 everything is exploratory.
- **Technical confounds are not modelled.** Fragment length, read extension, bigWig smoothing and
  anchors chosen from one condition can create width or shift changes. Fill in `bigwig_params` and
  `fragment_length` and use condition-independent anchors.
- **Background windows.** A neighbouring peak inside the outer flank inflates the background. Widen
  `flank` or shrink `background.width`.
- **Not implemented yet:** count models from BAM, segmentation, IWTomics, equivalence tests.
