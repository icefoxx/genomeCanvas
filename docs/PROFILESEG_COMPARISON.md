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

## 7. Ideas from FunChIP: shape classes and anchor alignment (evaluated)

FunChIP (Parodi et al. 2017, *Bioinformatics*, doi:10.1093/bioinformatics/btx201; Bioconductor until 3.18) clusters ChIP-seq
peaks by shape with a shift-only k-mean alignment of B-spline-smoothed profiles; it has no test
between conditions. We borrowed the clustering/alignment, always on the **pooled** profile of each
region (labels unused). Script: `benchmarks/profileseg_shapeclass_tuning.py` (tuning seeds
1 300 000+), validations with `benchmarks/profileseg_region_validation.py --variant`.

| variant | per-region calibration | per-region power (shift100 / shift50 / broaden) | aggregate centroid | decision |
|---|---|---|---|---|
| V0 current | max false share 0.0017, max FDP 0.054 | 0.80 / 0.20 / 0.33 (kernel 100, shift 50, kernel broaden) | null false calls up to 0.10 in four cells | default |
| V1 null stratified by shape class | false share 0, FDP 0.051 | lower: 0.75 / 0.15 / 0.27; frozen validation (seeds 1 800 000): calibrated, but 2 vs 2 real-template `shift` power 0.61 → 0.24 and an extra V4 failure | unchanged | **rejected**; low-level `strata` arguments kept for reproduction only |
| V2 anchors aligned | unchanged | unchanged | null false calls 0; power up in some cells (templates 3v3 0.6 → 0.8), down in one (templates 2v2 0.3 → 0.1) | **option** `shape.align_anchors` |
| kernel + derivative term (FunChIP d1) | unchanged | 0.38 (vs 0.80) | — | rejected |

Frozen validation of the alignment option (seeds 1 700 000, source hash `7e2d6be7c7eaf849`): V1–V3
pass in all six cells; V4 fails in the same cell as the default (real templates, 2 vs 2).

Real data (yeast MNase, GSE280671): with alignment the WT 2+2 nulls still give 0 calls in all three
splits, and 87 more ORFs become testable (1,555 → 1,642; double mutant 1,421 → 1,503 changed).
Shape classes were biologically sensible (a canonical class with maxima at +65 and +225 bp, i.e.
+1/+2 nucleosomes 160 bp apart, matching the measured repeat length; and an upstream-dominated class)
but both were affected alike (93 % and 91 %), so classes added no insight here. A class-wide median
of all regions' shifts was −0.4 bp yet "significant" (p ≈ 0.001) because most regions do not move;
class summaries now report direction balance and size among called regions instead, and are an
option (`shape.classes`).

Net: the pooled null already copes with heterogeneous shapes (a useful negative result), alignment
is a modest option for imprecise anchors and aggregate statements, and classes are descriptive.

## 8. In-silico insertion on real replicates, a 3 vs 3 real null, and conventional practice

Script: `benchmarks/profileseg_insertion.py`. Base experiments are splits of real replicates of one
condition: yeast MNase (GSE280671; six genotypes with four replicates each → 18 splits of 2 vs 2) and
human PRO-seq (GSE233930; DMSO controls of three degron lines, two replicates each → four 3 vs 3
splits with one replicate of each line per group, so line differences cancel — this is also the
3 vs 3 real null). Known changes were inserted into a random 20 % of the shape-testable regions in
the second group only (two random insertion sets per split); amplitude-only insertion is a negative
control for shape tests.

### 8.1 Results (BH 5 %; power / FDP)

Yeast MNase, 2 vs 2 (all 18 splits):

| method | null: share called / FWER | amplitude-only | shift 20 bp | shift 50 bp | shift 100 bp | broaden ×1.3 | broaden ×1.5 |
|---|---|---|---|---|---|---|---|
| `shape_kernel` | 0.017 / 0.17 | 0.017 | 0.14 / 0.12 | 0.88 / 0.12 | 0.99 / 0.12 | 0.99 / 0.13 | 0.99 / 0.13 |
| `shift` | 0 / 0 | 0 | **1.00 / 0.05** | 0.88 / 0.05 | 0.34 / 0.05 | 0.20 / 0.05 | 0.21 / 0.05 |
| v2 centroid (moderated t) | 0.0001 / 0.06 | 0.0001 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| Welch t on centroid | 0 / 0 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| global Wilcoxon over regions | **significant in 78 % of null splits** | 0.78 | 0.78 | 0.75 | 0.72 | 0.72 | 0.67 |

Yeast, excluding the three batch-aligned splits (see 8.2; post hoc, shown alongside):

| method | null | amplitude-only | shift 20 | shift 50 | shift 100 | broaden ×1.3 | broaden ×1.5 |
|---|---|---|---|---|---|---|---|
| `shape_kernel` | 0 / 0 | 0 | 0.06 / 0.03 | 0.85 / 0.04 | 0.99 / 0.04 | 0.99 / 0.04 | 0.99 / 0.04 |
| `shift` | 0 / 0 | 0 | 1.00 / 0.04 | 0.89 / 0.04 | 0.37 / 0.04 | 0.19 / 0.03 | 0.21 / 0.04 |

Human PRO-seq, 3 vs 3 (four line-balanced splits):

| method | null | amplitude-only | shift 20 | shift 50 | shift 100 | broaden ×1.3 | broaden ×1.5 |
|---|---|---|---|---|---|---|---|
| `shape_kernel` | 0 / 0 | 0 | 0.28 / 0.009 | 0.95 / 0.015 | 1.00 / 0.015 | 0.07 / 0.007 | 0.33 / 0.011 |
| `shift` | 0 / 0 | 0 | **0.97 / 0.006** | 1.00 / 0.006 | 0.41 / 0.006 | 0.24 / 0.004 | 0.52 / 0.004 |
| v2 centroid (moderated t) | 0 / **0.25** | 0 | 0.00 / 0.25 | 0.00 / 0.25 | 0.43 / 0.008 | 0.00 / 0.13 | 0.00 / 0.10 |
| Welch t on centroid | 0 / 0 | 0 | 0.00 | 0.02 | 0.32 | 0.00 | 0.00 |
| global Wilcoxon over regions | **significant in 100 % of null splits** | 1.00 | 1.00 | 1.00 | 1.00 | 0.75 | 1.00 |
| pausing index, Welch | 0 / 0 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| pausing index, moderated t | 0.0001 / **0.75** | 0.0005 | 0.00 / 0.63 | 0.04 / 0.07 | 0.27 / 0.04 | 0.00 / 0.44 | 0.00 / 0.28 |

### 8.2 Findings

- **Calibration on real noise.** `shift` had no false calls in any real null split (yeast 18, PRO-seq 4)
  and FDP ≈ 0.04–0.05 (yeast) and < 0.01 (PRO-seq). `shape_kernel` had none in PRO-seq and in 15 of
  18 yeast splits.
- **The three yeast splits with kernel calls are all "replicates 1+2 vs 3+4"** (double mutant 7.5 %,
  rtf1Δ 4.9 %, isw1Δ rtf1Δ 17.7 % of regions; all other splits 0). An independent check confirms
  batch-like pairing: in those genotypes the meta-profile distance across the pairs is 1.5–2.7× the
  distance within pairs, versus ~1.0 in WT and rtf1-8-9A (no calls). The kernel test detects real
  differences between replicate pairs when the split coincides with them; the directional test did
  not respond. Designs should balance processing batches across conditions; batch-aware label
  assignments (swaps only within batch) are a planned extension.
- **Complementary tests.** On sharp or periodic profiles the projection test is best for small
  shifts (20 bp: 0.97–1.00) but loses power for large ones (100 bp: 0.34–0.41: the projection is a
  local, linear test and a 100 bp move of a 160 bp nucleosome array is partly aliased), where the
  kernel test reaches 0.99–1.00. Read both; a region called by either is a shape change.
- **Conventional practice fails in characteristic ways.** Regions-as-units Wilcoxon calls 78–100 %
  of real null splits significant (pseudo-replication). Window-moment centroids have no power on
  real profiles (MNase arrays) and the moderated version is liberal on PRO-seq (FWER 0.25). The
  pausing index mixes position and amplitude: almost no power for shifts and FWER 0.75 with
  moderated variances.
- Amplitude-only insertion never produced shape calls beyond the null level.

## 9. ProfileSeg vs MMDiff2 on the same real BAMs

Script: `benchmarks/profileseg_realbam.py`. Yeast MNase-seq (GSE280671), first 3 M read pairs per
sample from ENA, bowtie2 to sacCer3, MAPQ ≥ 10, proper pairs (alignment rates 89–95 %). ProfileSeg
used fragment-midpoint counts (10 bp bins, ATG ± 1.5 kb) from these BAMs; MMDiff2 1.30.0 read the
same BAMs. MMDiff2's paired-end mode fails in this version (`getReads` references an undefined
`RawCounts.p`), so it analysed first mates in single-end mode with its own fragment-centre
estimation. A fixed random set of 300 genes (signal in every sample); ProfileSeg's shape-testable
subset of it (134–141 genes) is used for the like-for-like comparison.

| comparison (2 vs 2) | ProfileSeg `shape_kernel` | MMDiff2, BH < 0.05 (same genes) | MMDiff2 raw p < 0.05 |
|---|---|---|---|
| real null, WT 1+2 vs 3+4 | **0 / 134** | **17** | 0.19 |
| real null, WT 1+3 vs 2+4 | **0 / 134** | **7** | 0.11 |
| WT 1+2 vs isw1Δ chd1Δ 1+2 | **65 / 141** | 27 | 0.26 |

On all 300 genes MMDiff2 called 28 and 21 genes in the two nulls and 53 in the positive control.
This confirms on real data what the simulation and the source analysis indicated: MMDiff2's
p-values are anti-conservative in replicated designs, and its calls under a real change are only
about twice its false calls under the null. ProfileSeg made no null calls and detected 2.4× more
genes in the positive control. ProfileSeg on all 1,566–1,769 testable genes: 0 calls in both
nulls, 903 in the positive control (51 %, lower than the 91 % from full-depth GEO tracks because of
the 3 M-pair subsample).

Caveats: subsampled depth; MMDiff2 forced into single-end mode by its paired-end defect; one
dataset.

## 10. Runtime and memory

Scripts: `benchmarks/profileseg_resources.py` (components, simulation, 100 core bins, single thread,
fresh process per cell) and `benchmarks/profileseg_tool_resources.py` (whole tools on the same data).

### 10.1 ProfileSeg along the version timeline (seconds; current code)

| regions | design | v1 / v1.1 (global, 1000 bootstraps) | v2 (+ per-region amplitude) | v2.1 (+ kernel, template fit) | v2.2 (+ shift score) | v2.3 default (+ replicate diagnostics) |
|---|---|---|---|---|---|---|
| 1,000 | 2v2 | 0.3 | 0.3 | 0.4 | 0.4 | 0.6 |
| 1,000 | 3v3 | 0.4 | 0.4 | 0.6 | 0.6 | 1.0 |
| 5,000 | 2v2 | 1.2 | 1.3 | 1.9 | 1.9 | 2.8 |
| 5,000 | 3v3 | 1.7 | 1.8 | 2.7 | 2.7 | 4.4 |
| 5,000 | 4v4 | 2.2 | 2.2 | 3.6 | 3.7 | 7.9 |
| 20,000 | 3v3 | 6.8 | 6.9 | 10.8 | 11.1 | 17.8 |
| 20,000 | 4v4 | 9.3 | 9.5 | 16.4 | 16.8 | 36.4 |

Peak process memory: 0.4–0.6 GB at 1,000 regions, 0.8–1.5 GB at 5,000, 2.7–3.6 GB at 20,000
(8 samples). Options: shape classes and anchor alignment add ~5 s at 5,000 and ~20–24 s at 20,000
regions. Extracting profiles from bigWigs is not included (yeast: 24 tracks × 3,592 regions in ~1 min).

**Kernel feature map.** The kernel distance xᵀKx is computed as ‖Lᵀx‖² with K = LLᵀ (numerically
zero directions dropped), which gives identical p-values (unit test) at O(bins) instead of O(bins²)
per label assignment. Before this change the replicate diagnostics, which rerun the kernel test once
per left-out replicate, dominated:

| regions | design | kernel test, before → after | default pipeline, before → after |
|---|---|---|---|
| 5,000 | 3v3 | 4.5 s → 0.07 s | 32.7 s → 4.4 s |
| 20,000 | 4v4 | 77.8 s → 1.7 s | 658.6 s → 36.4 s |
| 5,000 | 6v6 (462 label assignments) | 441 s → 6.2 s | — |

The template fit is now processed in blocks of 2,000 regions (identical output; traced memory
680 → 227 MB at 6,000 regions, 4 vs 4).

### 10.2 Whole tools on the same simulated data (3 vs 3, single thread)

| tool | 300 regions: time / peak RSS | 1,000 regions: time / peak RSS |
|---|---|---|
| ProfileSeg (global with 1000 bootstraps + per region) | 3.3 s / 0.32 GB | 7.8 s / 0.40 GB |
| IWTomics (4 interval-wise tests, B = 1000) | 3.5 s / 0.64 GB | 4.0 s / 0.82 GB |
| edgeR QL diffSplice (2 count tables) | 1.2 s / 0.26 GB | 2.9 s / 0.17 GB |
| MMDiff2 (BAM → p-values) | 418 s / 7.8 GB | 912 s / 6.6 GB |

(ProfileSeg timings here predate the kernel feature map; they are upper bounds.) On the real yeast
BAMs (3 M read pairs per sample, 300 genes, 2 vs 2) MMDiff2 took ~15 s per comparison and ProfileSeg
0.2 s; MMDiff2's cost grows with reads per region. DANPOS was not run.

## 11. Closing the main gaps (batch handling, effect size, kernel normalisation, re-validation)

### 11.1 Batch-aware label assignments (`benchmarks/profileseg_batch_aware.py`)

The manifest `batch` column is now used. If condition and batch coincide (no assignment within
batches exists), the per-region shape tests are not performed (regions `not_tested`, with a note);
otherwise labels are permuted freely by default, or within batches with
`inference.batch_permutation: within`. Yeast MNase, batches = replicate pairs {1,2} / {3,4} (inferred
from the data; GEO records none), 2 vs 2, final code:

| split | batch mode | `shape_kernel` null share called | shift 50 bp: power / FDP | broaden ×1.3: power / FDP |
|---|---|---|---|---|
| confounded (1+2 vs 3+4) | off | 0.050 | 0.98 / **0.30** | 0.99 / **0.32** |
| confounded | default | not tested | — | — |
| balanced (1+3, 1+4 vs rest) | off = default | 0 | 0.68 / 0.017 | 0.97 / 0.017 |
| balanced | within | 0 | 0.93 / 0.060 | 0.99 / 0.059 |

### 11.2 `shift` effect size (`benchmarks/profileseg_shift_estimate_tuning.py`, tuning seeds 2 100 000+)

The pooled template of a region whose conditions differ is a blurred mixture. The template is now
rebuilt from samples re-aligned by their own estimated shifts (no labels; 3 refinements) with
7-bin smoothing (the score test keeps its validated 5-bin template). Worst-cell estimate/true ratio
0.58 → 0.73; real templates 2 vs 2 at 100 bp 0.70 → 0.86. In the final frozen validation V4 passes in
all cells (84.5 and 86.5 bp on real templates; 97.7–99.9 bp elsewhere). On PRO-seq the median
estimate is +18.8 / +19.3 bp (independent meta-profile ~+18 bp).

### 11.3 Kernel normalisation (found by the round-2 global validation)

At 4 vs 4 with 10 bp bins, `shape_kernel` called 1–3 % of regions under amplitude-only changes:
clipping the background-corrected profile at zero adds a noise floor whose relative weight grows as
the signal weakens, so an amplitude change looked like broadening. Dividing the signed profile by
its whole-window mass removed the bias but destroyed power on broad profiles (round-5 per-region
validation: real templates 2 vs 2, 100 bp shift power 0.59 → 0.02). Variants were compared on tuning
seeds (4 100 000+); the signed profile divided by its mass inside the condition-blind support
removes the bias (no calls under amplitude-only changes in any cell) and keeps power (real templates
2 vs 2: 0.62; parametric 2 vs 2: 0.73). This is the final implementation.

The exploratory per-region `log2_scale` test was removed: under translations its null FWER was
0.15–0.50.

### 11.4 Final frozen validations

Per-region (round 6, seeds 3 500 000+, source hash `d79dd8e1ec6c3d4f`): V1–V4 pass in all six cells.

| source | design | V2 false share (kernel / shift) | V3 FDP (kernel / shift) | V4 estimate (100 bp) | kernel power 100 / 50 / broaden | shift power 100 / 50 |
|---|---|---|---|---|---|---|
| parametric | 2v2 | 0 / 0 | 0.061 / 0.047 | 97.8 | 0.80 / 0.04 / 0.08 | 0.56 / 0.01 |
| parametric | 3v3 | 0 / 0 | 0.040 / 0.049 | 99.9 | 0.91 / 0.21 / 0.34 | 0.87 / 0.16 |
| semi-real | 2v2 | 0 / 0 | 0.046 / 0.056 | 97.7 | 0.72 / 0.04 / 0.24 | 0.86 / 0.16 |
| semi-real | 3v3 | 0 / 0 | 0.045 / 0.045 | 99.0 | 0.90 / 0.21 / 0.65 | 0.97 / 0.47 |
| real templates | 2v2 | 0 / 0 | 0.027 / 0.034 | 84.5 | 0.52 / 0.00 / 0.01 | 0.58 / 0.02 |
| real templates | 3v3 | 0 / 0 | 0.048 / 0.042 | 86.5 | 0.84 / 0.06 / 0.23 | 0.83 / 0.25 |

Global (`benchmarks/profileseg_validation.py`, round 4, seeds 390 000 / 395 000, source hash
`bee1aedc…`; criteria revised before round 2 as stated in the script: exact-null false calls,
Bonferroni across the two primary metrics):

| design | A1 | A2 | A3 | A4 | per-region shape tests (kernel / shift, max FDR R1 · R2) | per-region amplitude |
|---|---|---|---|---|---|---|
| V1 noise | 0.038 | 0.125 | **0.80** | **0.50**¹ | 0.040 · 0 / 0.025 · 0 | pass |
| V2a bimodal | 0.050 | **0.225** | **0.75** | **0.125** | 0.049 · 0 / 0.047 · 0.025 | `log2_auc` R2 0.15 |
| V2b skewed | 0.025 | 0.10 | 0.85 | 0.075 | 0.033 · 0 / 0.044 · 0 | pass (≤ 0.10) |
| V3a 3 vs 2 | 0.050 | 0.125 | 0.875 | 0.075 | 0.008 · 0 / 0.020 · 0 | absolute R2 0.20; relative R1 0.10, R2 0.15 |
| V3b 4 vs 4, 10 bp | 0.025 | 0.10 | 0.875 | 0.05 | 0.004 · 0.075 / 0.023 · 0.05 | `mad_width` R2 0.125 |
| V4 2-cell templates | 0.025 | 0.05 | **0.825** | 0.05 | 0.051 · 0 / 0.044 · 0.10 | absolute R2 0.275, R1 0.106 |

¹ In `global_with_site` no metric is null on average; A4 is evaluated per run on metrics whose truth
falls below the exact-null threshold in that run, so it rests on few runs (a separate run of that
scenario counted 9 eligible runs and gave 0). The criterion is fragile in such scenarios.

Reading: **the per-region shape tests pass every criterion in every design.** The remaining failures
are in (i) aggregate width metrics (`sd_width`/`iqr_width`: A2/A3 under amplitude changes — the known
re-weighting confound, flagged in reports) and (ii) per-region *absolute* amplitude in some designs
(`log2_signal`, `log2_auc`; relative amplitude in the 3 vs 2 design). These are documented
limitations; per-region absolute amplitude should be read with the diagnostics, and count-based
tools remain the reference for amplitude when BAMs are available.

### 11.5 Real data with the final code

| check | result |
|---|---|
| yeast WT 2+2 nulls (3) | 0 calls (both tests) |
| yeast insertion, batch-aligned splits excluded | `shape_kernel` FDP 0.032–0.036, power 0.73 (50 bp) / 0.97 (100 bp) / 0.96–0.98 (broaden); `shift` FDP 0.034–0.043, power 1.00 (20 bp) / 0.89 (50 bp) |
| PRO-seq 3 vs 3 line-balanced nulls (4) | 0 calls; insertion FDP ≤ 0.018 |
| yeast, same real BAMs vs MMDiff2 (134 genes) | nulls: ProfileSeg 0 / 0, MMDiff2 (300 genes) 33 / 20; positive: ProfileSeg 90 of 141 |
| yeast isw1Δ chd1Δ vs WT (4 vs 4) | 1,264 / 1,555 genes change shape; `shift` 510 |
| PRO-seq NELFB 2 h / 48 h (2 vs 2) | `shift` 6,350 / 2,834 genes, 98 % downstream, +18.8 / +19.3 bp; `shape_kernel` 7,990 / 5,166 (48 h was 0 with the clipped kernel) |

### 11.6 TF CUT&RUN: single-cell-line 3 vs 3 nulls and a second real-BAM comparison (GSE291134)

`benchmarks/profileseg_cutrun_null.py`, `benchmarks/profileseg_realbam_cutrun.py`. EWS::FLI1 CUT&RUN,
untreated, six replicates in each of two cell lines (A673; EF-dTAG clone A2.2). Anchors: 5,000
condition-blind summits of the pooled signal.

| check | result |
|---|---|
| A673, 10 splits 3 vs 3 (bigWigs) | `shape_kernel` and `shift`: 0 calls in every split; aggregate report: no calls |
| A2.2, 10 splits 3 vs 3 | 0 calls except one split (`shape_kernel`, 0.14 % of regions); aggregate: no calls |
| insertion on these replicates | power ≤ 0.12 even for 100 bp shifts |
| real BAMs (first 4 M pairs, chr1/chr2), 2 nulls + read-level insertion (+60 bp, 60 of 295 anchors) | ProfileSeg: 0 calls in both nulls, 0 of 60 inserted shifts detected; MMDiff2: stops with `KernelMatrix[posB, posB]: only 0's may be mixed with negative subscripts` (peaks with no reads in a sample) |

Fragments per anchor core (±600 bp) in these BAMs: median 4–20 per sample, versus 412 per gene in the
yeast MNase BAMs at a similar read depth. The FLI1 sites spread the reads over thousands of loci, so
this dataset carries too little per-site shape information at this depth. It confirms calibration
(single-line 3 vs 3 nulls, both per-region and aggregate) but leaves the second real-BAM comparison
inconclusive. Practical guidance: check the number of fragments per region per sample before
per-region shape testing — hundreds worked (yeast), around ten did not.
