# ProfileSeg compared with existing methods

Two parts: a methodological comparison (what each method estimates and how it does inference),
then an empirical comparison on simulations with known truth, with every tool run as published
(IWTomics 1.26.0, MMDiff2 1.30.0, edgeR on Bioconductor 3.18 / R 4.3.2). Scripts:
`benchmarks/profileseg_comparison.py` (pointwise baselines, edgeR) and
`benchmarks/profileseg_external.py` (IWTomics, MMDiff2).

## 1. Methodological comparison

| | ProfileSeg | IWTomics | MMDiff2 | pointwise t + BH | edgeR / DiffBind / csaw |
|---|---|---|---|---|---|
| question | does amplitude **or** shape change: typical region (global) and which regions (per region) | where along the aligned domain do two groups of curves differ | which peaks differ in read-position distribution | at which bins do mean profiles differ | which regions/windows differ in read count |
| input | coverage tracks (bigWig) or computeMatrix | aligned feature matrices for two region sets | BAMs + peaks | aggregate profiles | BAM counts |
| unit of inference | **biological replicate** | the curve: a **region** in its intended use (two region sets) | peak; replicates used to build the null | whatever the curve is | replicate |
| replicate variability | Welch on replicate summaries; per region, moderated t (variance shrinkage over regions) | not modelled when curves are regions; with replicates as curves the permutation space is tiny (3 vs 3: 20 assignments, p ≥ 0.1) | empirical null from within-condition replicate MMDs, smoothed over coverage (locfit) | replicates as observations; no sharing across bins | negative binomial GLM with dispersion shrinkage |
| amplitude vs shape | separated by construction: background-corrected, mass-normalised shape metrics; amplitude metrics reported apart | not separated; the user has to normalise curves | intended to be shape-only, but MMD depends on coverage | not separated unless curves are area-normalised | amplitude only |
| estimand | interpretable: log2 FC, centroid shift (bp), spread change (bp), asymmetry | a significance map, no effect size | one distance per peak, no direction | per-bin difference | log2 FC |
| multiplicity | Holm over secondary global metrics; BH over regions | interval-wise FWER over the domain | BH over peaks | BH over bins | BH over regions |
| smallest p, 3 vs 3 | Welch: continuous; exact permutation floor 0.1 reported | 0.1 when curves are replicates | continuous (parametric locfit null) | continuous | continuous |
| strength | separates confounds; effect sizes in bp; self-benchmark on the user's own shapes | localises *where* a difference sits, at all scales | uses the full read distribution; sensitive to complex shape changes | simple | mature count models, spike-in support |
| weakness | per-region shape has low power; summary metrics miss shape changes that leave centroid and spread unchanged (e.g. a symmetric change of kurtosis) | pseudo-replication when regions are curves; no power when replicates are curves | coverage-dependent distance and null; no effect size or direction; needs BAMs | confounded by amplitude and background; no effect sizes | blind to shape |

The key methodological point is the **unit of inference**. A WT-vs-KO experiment has a few replicates per
condition and hundreds to thousands of regions measured in every sample. Methods that treat regions
(IWTomics in its intended use) or reads (the plain MMD bootstrap, which the MMDiff authors found
rejects nearly every replicate-vs-replicate comparison) as independent units test whether *these
particular samples* differ, not whether the *conditions* differ. ProfileSeg and MMDiff2 both use
replicate variability, ProfileSeg directly (replicate summaries, moderated t) and MMDiff2 through
an empirical null.

## 2. Empirical comparison

Simulations: 3 vs 3 replicates, 300 regions, negative-binomial counts with replicate and library
noise (`profileseg.simulate`). Every method sees the same data: IWTomics gets the core ±1000 bp
in 40 bp bins (per region, or per replicate as the mean over regions; raw or area-normalised);
MMDiff2 gets single-end BAMs (36 bp reads, 200 bp fragments) simulated from the same counts.
`shape changes = no` means any shape call is a false positive.

### 2.1 Global: share of runs reporting a profile/shape difference (10 runs per scenario)

| scenario | shape changes | IWTomics regions, raw | IWTomics regions, area | IWTomics replicates (raw or area) | ProfileSeg |
|---|---|---|---|---|---|
| null | no | 0.00 | 0.00 | 0.00 | 0.00 |
| amplitude (unequal across regions) | no | **1.00** | **1.00** | 0.00 | 0.00 |
| global_scale | no | **1.00** | 0.10 | 0.00 | 0.00 |
| site_amplitude (20 % of regions ×2.5) | no | **1.00** | **1.00** | 0.00 | 0.10 |
| shift | yes | 1.00 | 1.00 | 0.00 | 1.00 |
| broaden | yes | 1.00 | 1.00 | 0.00 | 1.00 |
| opposite_shift | yes | 1.00 | 1.00 | 0.00 | 1.00 |

With replicates as curves the smallest attainable IWTomics adjusted p is 0.1 (3 vs 3: 20 assignments).
With regions as curves IWTomics answers "do these two sets of curves differ", which any
amplitude change makes true. Area normalisation removes a uniform rescaling but not
region-to-region amplitude changes, because those re-weight the mean curve.

The earlier comparison (`benchmarks/profileseg_comparison.py`, 40 runs) gives the same picture
for simple pointwise baselines: pointwise t-tests with regions as units call the amplitude-only
scenario in 100 % of runs (replicates as units: also 100 %, because raw mean profiles carry the
amplitude change), and ProfileSeg calls it in 5 %. A fragment-length change (a real but subtle shape
change) is detected in 30 % of runs by ProfileSeg and 5–25 % by the pointwise methods.

### 2.2 Per region: share of regions called (BH < 0.05 unless stated)

| scenario | shape changes | MMDiff2 (BH < 0.05) | MMDiff2 (own padj < 0.1) | MMDiff2 calibrated* | ProfileSeg centroid | ProfileSeg mad_width | ProfileSeg either |
|---|---|---|---|---|---|---|---|
| null | no | 0.34 | 0.38 | 0.05 | 0.00 | 0.00 | 0.00 |
| amplitude | no | 0.81 | 0.84 | 0.44 | 0.00 | 0.00 | 0.00 |
| global_scale | no | 0.34 | 0.38 | 0.05 | 0.00 | 0.00 | 0.00 |
| site_amplitude | no | 0.52 | 0.56 | — | 0.00 | 0.00 | 0.00 |
| shift (100 bp, every region) | yes | 0.99 | 0.99 | 0.83 | 0.40 | 0.00 | 0.40 |
| broaden (every region) | yes | 0.68 | 0.72 | 0.23 | 0.00 | 0.02 | 0.02 |
| opposite_shift (every region) | yes | 0.99 | 0.99 | 0.78† | 0.39 | 0.00 | 0.39 |

\* Raw p-value threshold chosen post hoc so that 5 % of regions are called in the null runs
(p ≤ 2.8e-10); this is an oracle calibration that a user cannot do, used only to compare power at
equal false-call rate. † From one run; site_amplitude was not available when the calibration was
computed (the script now reports this column directly).

**Why MMDiff2 is anti-conservative here.** The simulation is exchangeable under the null (every
sample's noise is independent), so the 34 % null call rate comes from the p-value construction in
MMDiff2 1.30.0 (`compute.pvalues` / `detpvals`): for each peak, the 9 between-group pairwise
distances (3 × 3, sharing samples, hence dependent) each get a p-value from a locfit null, and the
9 p-values are combined with Fisher's method, which assumes independence; the dispersion is fitted
on positive residuals only; and the t reference uses the number of locfit points as degrees of
freedom. Under the null 24 % of raw p-values were below 0.001. This is a finding on simulated BAMs
and has not been checked on real data.

**MMDiff2 also responds to amplitude.** Its distance compares read-position distributions without
removing background, so a change in peak-to-background ratio changes the distribution; even after
calibration 44 % of regions are called in the amplitude-only scenario. A MMDiff2 call therefore
means "the read distribution differs", not "the peak shape changed".

**ProfileSeg per region** has essentially no false calls in these scenarios but limited power: a
100 bp shift is found in 40 % of regions (MMDiff2 calibrated: 83 %), and spread changes are almost
never found per region (2 % vs 23 %). Spread changes are reliably detected only at the aggregate
level. MMDiff2 uses every read position and has more information per region; ProfileSeg reduces
each region to a few summaries.

**MMDiff2 implementation note.** With 200 or more peaks, `compDists` samples
`numPeaksForSigma = 100` peaks from those strictly inside the count interquartile range to set the
kernel bandwidth; when fewer than 100 qualify it stops with "cannot take a sample larger than the
population". The benchmark uses 300 peaks to avoid this; with smaller peak sets the bandwidth has to
be supplied by hand.

### 2.3 Per-region amplitude (40 runs, `benchmarks/profileseg_comparison.py`)

| scenario | method | FDR | power |
|---|---|---|---|
| site_amplitude | Welch per region | 0.016 | 0.15 |
| site_amplitude | moderated t (ProfileSeg) | 0.035 | 0.99 |
| site_amplitude | edgeR QL + TMM on window counts | 0.063 | 0.99 |
| global_with_site | Welch per region | 0.010 | 0.10 |
| global_with_site | moderated t (ProfileSeg) | 0.046 | 0.88 |
| global_with_site | edgeR QL + TMM | 0.069 | 0.88 |
| amplitude (no site-specific change) | edgeR QL + TMM | 0.225 (0.4 discoveries/run) | — |

For amplitude, ProfileSeg's per-region test and edgeR have the same power; edgeR on window counts
(which include background) is slightly anti-conservative when IP efficiency changes, which the
background-corrected signal avoids. With BAMs and a count model available, DiffBind/csaw remain
the reference for amplitude.

### 2.4 Caveats on fairness

- The data come from ProfileSeg's own simulator. The scenarios are chosen to separate amplitude
  from shape, which is the question ProfileSeg was built for; they favour it by construction.
- IWTomics is used outside its intended design (two region sets) when applied to WT vs KO on the
  same regions; that mismatch is the point being tested, not a flaw of IWTomics.
- MMDiff2 gets simulated single-end BAMs whose reads are placed from binned counts; real BAMs have
  fragment-length and mappability structure this does not reproduce.
- 10 runs per scenario: run-level shares have a binomial SE of up to 0.16.
- ProfileSeg's summaries miss shape changes that leave centroid and spread unchanged (e.g. a
  symmetric change in kurtosis or a bimodal split); MMDiff2 and IWTomics would detect those.

### 2.5 When to use which

| question | tool |
|---|---|
| which peaks change in read count | DiffBind / csaw / edgeR (BAMs), or ProfileSeg `log2_signal_relative` (coverage tracks) |
| does the typical profile change amplitude or shape, with replicates, separated | ProfileSeg global |
| which regions shift, with an effect in bp | ProfileSeg per region (large shifts only), MMDiff2 for screening with its p-values treated as rankings |
| where along the window two **region sets** differ | IWTomics |
| any read-distribution difference at a peak, including complex shapes | MMDiff2, with external calibration (e.g. permutation of sample labels) and an amplitude check |

## 3. Improving per-region shape power: candidate methods

`benchmarks/profileseg_region_methods.py`, 20 runs per scenario, 3 vs 3, 300 regions; parametric
shapes and real (GSE112560 TCF3) templates. `mix_*`: 20 % of regions change (the realistic case).
Cells: share of regions called (no-change scenarios) or power / FDP (change scenarios), BH < 0.05.

| method | false calls (null / amplitude / site) | shift: param / real | mix_shift: param / real | broaden: param / real | mix_broaden: param / real | effect estimate (true shift 100 bp) |
|---|---|---|---|---|---|---|
| A0 current centroid / MAD, full window | ≤ 0.001 | 0.58 / 0.00 | 0.11 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 100 / 62 bp |
| A support-restricted centroid / MAD | ≤ 0.001 | 1.00 / 0.04 | 0.99 / 0.00 | 0.64 / 0.00 | 0.56 / 0.00 | 67 / 40 bp (biased low) |
| B A + precision weights | ≤ 0.007 | 1.00 / 0.06 | 0.99 / 0.00 | 0.67 / 0.00 | 0.56 / 0.00 | as A |
| C1 edgeR QL diffSplice, all bins | **0.48 / 0.20** under amplitude / site (param) | 0.95 / 0.97 | 0.77 / 0.88 | 0.55 / 0.76 | 0.18 / 0.45 | none |
| C1 edgeR QL diffSplice, support bins | 0.000 | 0.95 / 0.93 | 0.84 / 0.82 | 0.43 / 0.50 | 0.16 / 0.17 | none (omnibus) |
| C2 template fit (shift / log scale) | ≤ 0.001 | 0.93 / 0.46 | 0.79 / 0.13 | 0.19 / 0.00 | 0.07 / 0.00 | 88 / 84 bp |
| C2w C2 + precision weights | ≤ 0.002 | 0.96 / 0.66 | 0.88 / 0.37 | 0.28 / 0.01 | 0.04 / 0.00 | 89 / 82 bp |
| D groups of 20 regions (C2 means) | shift ≤ 0.01; **scale 0.14 under amplitude** | 1.00 / 0.98 (groups) | 1.00 / 0.90 | 1.00 / 0.07 | 1.00 / 0.13 | 87 / 81 bp |
| E kernel statistic, pooled permutation null | 0.000 | 0.99 / 0.96 | 0.85 / 0.72 | 0.41 / 0.36 | 0.07 / 0.04 | none (omnibus) |

Findings:
- The main loss in the current method is the estimator, not the data: on real templates a 100 bp
  shift goes from 0 % (full-window centroid) to 72–96 % with tests that use the whole shape (E, C1).
  Moment summaries of broad, irregular real profiles carry little of the shift information.
- Restricting moments to the condition-blind support (A) fixes parametric shapes but not real ones,
  and attenuates the effect estimate; precision weights (B) add nothing and slightly inflate false calls.
- Count GLMs on all bins (C1 all) confuse amplitude with shape (as MMDiff2 does): background bins
  change share whenever the peak changes height. Restricting to support bins removes this.
- The kernel test (E) is the best coverage-track method: no false calls, high power for shifts on
  both shape sets, no counts needed. It is omnibus (no direction or bp); pair it with C2 for effect sizes.
- Spread changes in 20 % of regions remain hard for every method (≤ 0.17 on real templates).
- Group tests (D) have high power but answer a different question, and the scale version is
  confounded by amplitude.

Caveats: the simulator draws bin counts independently given the rates, which is exactly edgeR's
model, so C1 is favoured; with real data C1 must use fragment-midpoint counts from BAMs (coverage
counts one fragment in several bins and would make it anti-conservative). E's pooled null assumes the
normalised statistic has the same null distribution across regions; both need a real-data null check
(replicate-vs-replicate with ≥ 6 samples of one condition). With 2 vs 2 E has only 3 label splits.

## 4. v2.1: tuning and frozen validation of the per-region shape tests

### What changed in `per_region_tests`

- `shape_kernel` (new, omnibus): Gaussian-kernel distance between condition-mean mass-normalised
  profiles over the full core window (σ = 60 bp), divided by within-condition spread plus s0,
  with a null pooled over regions and all non-observed replicate-label assignments. Works for any
  n1 vs n2 (2 vs 2: 2 null assignments per region). Calls: Changed / Not_significant.
- `shift` / `log2_scale` (new): per-sample template fit (weighted Gauss–Newton on the pooled,
  lightly smoothed profile); moderated t. `shift` is the effect size in bp.
- `centroid` / `mad_width`: unchanged full-window moments (support-restricted moments were
  miscalibrated on real noise: null FWER 0.2–0.4 in tuning).
- Shape tests now use every region with positive pooled mass (condition-blind).

### Tuning (`benchmarks/profileseg_region_tuning.py`, seeds 100 000+)

Sources: parametric, real templates, semi-real (`benchmarks/semireal.py`: real GSE112560 replicate
residuals added to a known truth) × designs 3v3 and 2v2. Three rounds:

1. The pre-specified rule lacked an FDP constraint; its pick had FDP 0.29 in one cell. The rule was
   amended (per-cell FDP ≤ 0.10, matching validation criterion V3) **before any validation run**.
2. Scaling the pooled null by a label-invariant spread did not help (rejected).
3. Root cause found: s0 was a quantile of W computed **per label assignment**. Assignments that mix
   conditions carry part of a true effect in W, so the null assignments got a larger s0 than the
   observed one, making the test anti-conservative once many regions change (in one dataset 62
   false calls out of 120 vs 8 with a common s0). s0 is now a quantile of the total spread, which
   is identical under every assignment. Selected: full window, σ = 60 bp, s0 = median; template
   smoothing 5 bins.

### Validation (`benchmarks/profileseg_region_validation.py`, seeds 300 000+, source hash `017027fc42d30416`)

End to end through `per_region_tests`; 20 runs per cell; criteria fixed before running:
V1 kernel null FWER ≤ 0.10; V2 false-call share ≤ 0.005 in no-shape scenarios (kernel and shift);
V3 mean FDP ≤ 0.10 in 20 % mixtures (kernel and shift); V4 shift estimate within ±20 % of 100 bp.

| source | design | V1 | V2 (kernel / shift) | V3 (kernel / shift) | V4 shift estimate | result | kernel power: shift 100 / 50 / broaden | shift power: 100 / 50 |
|---|---|---|---|---|---|---|---|---|
| parametric | 3v3 | 0.00 | 0.0013 / 0.0000 | 0.039 / 0.009 | 88.0 | pass | 0.91 / 0.21 / 0.48 | 0.78 / 0.01 |
| parametric | 2v2 | 0.00 | 0.0000 / 0.0003 | 0.039 / 0.019 | 79.1 | **fail V4** | 0.73 / 0.03 / 0.17 | 0.37 / 0.00 |
| templates | 3v3 | 0.00 | 0.0000 / 0.0000 | 0.044 / 0.017 | 81.9 | pass | 0.87 / 0.06 / 0.35 | 0.18 / 0.00 |
| templates | 2v2 | 0.00 | 0.0000 / 0.0010 | 0.038 / **0.200** | 72.4 | **fail V3, V4** | 0.54 / 0.00 / 0.03 | 0.00 / 0.00 |
| semi-real | 3v3 | 0.00 | 0.0029 / 0.0000 | 0.062 / 0.045 | 99.9 | pass | 0.78 / 0.05 / 0.23 | 0.83 / 0.09 |
| semi-real | 2v2 | 0.00 | 0.0000 / 0.0012 | 0.075 / 0.034 | 97.3 | pass | 0.60 / 0.02 / 0.17 | 0.61 / 0.00 |

Reading:
- `shape_kernel` met V1–V3 in all six cells, including 2 vs 2 and real replicate noise. On
  real-noise data its FDP (0.06–0.075) is above the nominal π0·α ≈ 0.04 but within the criterion.
- `shift` failed in two 2 vs 2 cells: its estimate is attenuated by 21–28 % (pre-registered
  tolerance 20 %), and with real templates at 2 vs 2 it has no power, so its few calls are mostly
  false (mean FDP 0.20 from occasional single false calls). These failures stand; `shift` should
  be reported as an effect size with that attenuation stated, and not used for calling at 2 vs 2
  on broad profiles.
- Power for 50 bp shifts and for broadening in 20 % of regions remains low for every method
  (kernel ≤ 0.21 and ≤ 0.48).
- The count GLM (C1) is not integrated: it needs fragment-midpoint counts from BAMs.

### v2.1 against MMDiff2 on the same datasets

ProfileSeg v2.1 re-run on the exact datasets of §2.2 (same seeds; 3 vs 3, 300 regions, 10 runs;
every region changes in the shape scenarios). Share of regions called:

| scenario | shape changes | MMDiff2 (BH < 0.05) | MMDiff2 calibrated* | ProfileSeg v2 (centroid or MAD) | **v2.1 `shape_kernel`** | v2.1 `shift` |
|---|---|---|---|---|---|---|
| null | no | 0.34 | 0.05 | 0.00 | **0.00** | 0.00 |
| amplitude | no | 0.81 | 0.43 | 0.00 | **0.00** | 0.00 |
| global_scale | no | 0.34 | 0.05 | 0.00 | **0.00** | 0.00 |
| site_amplitude | no | 0.52 | — | 0.00 | **0.00** | 0.00 |
| shift 100 bp | yes | 0.99 | 0.83 | 0.40 | **0.96** | 0.93 |
| broaden ×1.5 | yes | 0.68 | 0.23 | 0.02 | **0.67** | 0.00 |
| opposite_shift | yes | 0.99 | 0.78 | 0.39 | **0.95** | 0.93 |

\* oracle threshold giving 5 % called regions under the null (not available to a user).

`shape_kernel` reaches MMDiff2's uncalibrated sensitivity with no false calls, and exceeds MMDiff2
at an equal per-region false-call rate. Caveats: parametric shapes from ProfileSeg's own simulator
(the v2.1 tuning used other seeds of the same simulator family), MMDiff2 on simulated single-end
BAMs, and different error metrics (BH FDR vs a per-region false-call rate). A comparison on real BAMs
with known positive controls is still missing.

## 5. v2.2: small shifts

### Why 50 bp shifts are hard (`benchmarks/profileseg_small_shift.py`)

Per-sample shift estimates (projection on the template derivative) have a replicate SD of 12–19 bp.
The realised z of a 50 bp shift in the shifted regions has a median of only 2.0–3.1 (3 vs 3 and
2 vs 2), while surviving BH over 300 regions needs roughly z > 3.5. Most of the limit is therefore
information (replicates, depth, peak width), not the method. The remaining loss was the omnibus
kernel test spending many degrees of freedom on a one-dimensional change. Multi-scale kernels,
lower s0 and grouped weighted BH (IHW-style) did not help (weighted BH raised null FWER to 0.1–0.2).

### Change

`shift` calls now come from `shift_score_test`: the projection of each sample's normalised profile on
the derivative of the pooled template (the locally most powerful direction for a small
translation), with the pooled permutation null and a label-invariant s0. The estimate (bp) still
comes from the template fit.

The semi-real generator was also corrected: plain residuals from condition means are dependent
(with two samples exactly r and −r); it now uses Helmert contrasts (four independent real noise
vectors) and draws reused vectors from donor regions of similar strength. The semi-real rows of §4
are superseded.

### Validation round 2 (seeds 700 000+, source hash `224924527b2cc935`, same criteria V1–V4)

| source | design | V1 | V2 (kernel / shift) | V3 (kernel / shift) | V4 shift estimate | result | kernel power: 100 / 50 / broaden | **shift power: 100 / 50** (round 1) |
|---|---|---|---|---|---|---|---|---|
| parametric | 3v3 | 0.00 | 0.0008 / 0.0005 | 0.030 / 0.039 | 89.3 | pass | 0.92 / 0.17 / 0.45 | **0.87 / 0.16** (0.78 / 0.01) |
| parametric | 2v2 | 0.00 | 0 / 0 | 0.041 / 0.036 | 80.8 | pass | 0.72 / 0.04 / 0.15 | **0.53 / 0.00** (0.37 / 0.00) |
| templates | 3v3 | 0.00 | 0 / 0 | 0.044 / 0.042 | 84.3 | pass | 0.88 / 0.07 / 0.33 | **0.83 / 0.23** (0.18 / 0.00) |
| templates | 2v2 | 0.00 | 0 / 0 | 0.052 / 0.050 | 72.8 | **fail V4** | 0.68 / 0.00 / 0.02 | **0.61 / 0.04** (0.00 / 0.00) |
| semi-real | 3v3 | 0.00 | 0 / 0 | 0.041 / 0.048 | 95.9 | pass | 0.91 / 0.21 / 0.61 | **0.97 / 0.46** |
| semi-real | 2v2 | 0.00 | 0 / 0 | 0.055 / 0.066 | 95.5 | pass | 0.76 / 0.03 / 0.18 | **0.89 / 0.19** |

Calibration criteria V1–V3 pass in every cell for both tests. The only failure is the template-fit
effect size at 2 vs 2 on real templates (attenuated by 27 %, tolerance 20 %), as in round 1. Round-1
semi-real numbers are not comparable (old generator).

## 6. Real-data validation

Scripts: `benchmarks/profileseg_replicate_artefact.py` (simulation), and for real data the
analyses described below (bigWigs from GEO, read with pyBigWig only). Code state: v2.2 plus the
replicate-influence diagnostics.

### 6.1 One replicate with a sample-wide artefact (simulation)

A 30 bp shift or ×1.25 broadening of every region in one replicate did **not** increase false
calls or FDP in any source or design (every label assignment contains the deviant replicate, so the
pooled null absorbs it). The group-aware offset check flagged the artefact in 100 % of runs and named
the right replicate (3 vs 3) or the right pair (2 vs 2), with no false alarms in clean data. An
earlier version compared each replicate with all samples and raised false alarms in up to 60 % of
runs with genuine condition effects; it was replaced.

### 6.2 Yeast MNase-seq: real nulls and a known positive control (GSE280671)

Mononucleosome CPM bigWigs; 3,592 ORFs ≥ 1 kb, anchored at the ATG (oriented), ±1.5 kb, 10 bp bins.
Known biology: isw1Δ chd1Δ disrupts nucleosome spacing and phasing genome-wide; single mutants
have milder defects.

| comparison | `shape_kernel` calls | `shift` calls |
|---|---|---|
| null: WT 2+2 splits (3 splits) | **0, 0, 0** of 1,490 | **0, 0, 0** |
| null: double-mutant 2+2 splits (3 splits) | 0, 0, 122 (7.5 %) of 1,624 | 0, 0, 0 |
| WT (4) vs chd1Δ (2) | 789 / 1,383 | 482 |
| WT (4) vs isw1Δ (2) | 859 / 1,506 | 421 |
| WT (4) vs isw1Δ chd1Δ (4) | 1,421 / 1,555 (91 %) | 510 |
| WT (2) vs isw1Δ chd1Δ (2), three WT pairs | 92–95 % | 410–445 |

Independent check (meta-profile autocorrelation at the nucleosome repeat, a phasing measure not used
by ProfileSeg): WT 0.96–0.97, single mutants 0.89–0.94, double mutant 0.40–0.67. The ordering of
effects agrees with ProfileSeg and with the literature. The double-mutant replicates are themselves
heterogeneous (meta-profile distances 0.03–0.11 vs 0.016–0.03 among WT; phasing 0.67 / 0.54 / 0.40 /
0.61), and the one null split with calls (rep1+2 vs rep3+4) is the split aligned with that gradient;
the other two splits give none. We read the 122 calls as a genuine replicate difference, not a
method error, but this cannot be proven from the data.

The leave-one-out flag fired on an unremarkable WT replicate in 2 vs 2 positives (dropping one
replicate from 2 vs 2 mainly costs power); it is now restricted to designs with ≥ 3 replicates per
group.

### 6.3 Human PRO-seq: NELF degradation (GSE233930), 2 vs 2

RefSeq Select protein-coding TSSs (main chromosomes, genes ≥ 2 kb), strand-matched sense signal,
±1 kb, 10 bp bins. Known biology: loss of NELF moves promoter-proximal Pol II downstream.

| | independent meta-profile (top 2,000 promoters) | `shift` calls (direction) | median shift estimate | `shape_kernel` calls |
|---|---|---|---|---|
| DMSO | pause peak 35 bp; centroid (0–300 bp) 58–64 bp | — | — | — |
| NELFB degraded 2 h | peak 45 bp; centroid 80–81 bp | 6,350 / 13,797 (**98 % downstream**) | **+17.5 bp** | 4,223 |
| NELFB degraded 48 h | peak 45 bp; centroid 74–81 bp | 2,834 / 13,923 (**98 % downstream**) | **+18.4 bp** | **0** |

The direction (98 % downstream; false calls would be symmetric) and the size (+17–18 bp vs ~18 bp
in the independent meta-profile) agree with the known biology. At 48 h the omnibus kernel test found
nothing while the directional test found 2,834 genes: an 18 bp shift is small relative to the 60 bp
kernel, as predicted in §5 — the two tests are complementary. The aggregate centroid (+18.8 bp)
was "no clear change" at 2 vs 2 (Welch interval with 1–2 df), whereas the per-region tests, which
borrow information across regions, detected it: aggregate statements need ≥ 3 replicates.

### 6.4 What remains

- No real null with ≥ 6 replicates of one condition (3 vs 3 nulls) was available as per-replicate
  tracks; the yeast nulls are 2 vs 2 and the splits are not independent.
- Comparison with MMDiff2 on these real BAMs was not possible (only bigWigs are deposited).
- Positive controls show direction and size agree with biology; they cannot give per-region FDR.
