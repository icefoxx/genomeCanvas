# ProfileSeg: implementation and benchmark specification (v2)

> Revision of `ProfileSeg_Agent_Implementation_Plan.md`. Main changes: background/SNR confounding
> of shape metrics, fragment-length and smoothing confounds, anchor-selection bias, orientation
> rules, condition-independent filtering, a concrete inference plan for 2–3 replicates, clarified
> units of inference for IWTomics/DiffSegR, semi-synthetic benchmarks, and a staged scope.
> Statements about IWTomics and DiffSegR internals must be checked against their documentation at
> Gate 1; they are marked *(verify)*. v1 is implemented as `genomecanvas.profileseg`
> (`gcanvas-profileseg`); see docs/PROFILESEG.md and the "v1 implementation findings" in Module 2.

## Mission
Develop a reproducible **prototype**, not a claimed validated method, for comparing aligned genomic
signal profiles from ChIP-seq, CUT&RUN and ATAC-seq between conditions. Combine:
- amplitude summaries;
- shape/transport metrics computed on **background-corrected** signal;
- (later) change-point segmentation;
- statistically defensible, **replicate-level** inference.

Do not claim that DiffSegR handles aligned bigWig profiles, or that FPOP segmentation produces
p-values.

**Staging.**
- **v1:** extraction, QC, amplitude, shape metrics, replicate-level effect sizes with uncertainty,
  semi-synthetic benchmark, figures A/B/C/E/F.
- **v2:** segmentation (Module 3), local inference, count-based (BAM) models, IWTomics comparison.

Each stage must pass its gate before the next starts.

## First actions: environment, data, metadata (Gate 1)
1. Read the local skills/instructions directory if one exists (e.g. `~/scripts/skills/`). Report
   exact filenames and commands; do not assume interfaces. This step is environment-specific.
2. Search the project for bigWigs, BAMs, peak BEDs, computeMatrix files, sample sheets, replicate
   metadata and reports. Produce a discovery report before any expensive job.
3. Write an input manifest with these columns (blank where unknown; flag missing values rather
   than guessing):
   `sample_id, condition, biological_replicate, assay, genome_build, signal_file, signal_type,
   normalization, region_set, batch, pairing_id, bigwig_params, fragment_length`.
   - `bigwig_params` records how coverage was made: bin size, `--extendReads`, `--smoothLength`,
     offsets/Tn5 shift, normalisation. Signal shape depends on these.
   - `fragment_length` is the median insert size (from BAM or library QC).
4. Confirm the genome build and coordinate conventions; never mix builds. Handle missing, zero
   and negative signal separately.
5. Check versions and primary documentation for any comparator actually used (IWTomics, DiffSegR,
   FPOP implementations, DESeq2/edgeR, limma). Cite them and note deviations.

## Scientific questions: distinct estimands
- **A. Amplitude:** mean coverage, AUC and replicate-aware contrast on the **same fixed regions**.
- **B. Shape:** differences after background correction and area normalisation; transport and
  distributional statistics.
- **C. Shift/width:** signed centroid displacement (oriented anchors only) and change in spatial
  spread (SD and robust/IQR width).
- **D. Local effects:** intervals where the contrast differs, with controlled error (v2).
- **E. Heterogeneity:** the distribution of region-level changes, not only group means.

These are not interchangeable:
- Area normalisation removes global multiplicative changes.
- A large W1 does not identify a peak shift.
- Before looking at results, pre-specify **one primary metric per question** (e.g. ΔAUC for A,
  Δcentroid for C). The other metrics are secondary, so the number of metrics is itself a
  multiplicity issue.

## Input and extraction
- **Inputs:** (1) replicate bigWigs + BED regions; (2) precomputed region × bin × replicate
  matrices (deepTools `computeMatrix` or `.npz`); (3, v2) BAMs for raw counts.
- **Anchors must not be chosen from one condition.** Summits called in WT and then applied to KO
  create apparent KO shifts through regression to the mean. Use condition-independent anchors:
  motif or TSS positions, or summits/regions called on **pooled** data (union or consensus).
  Record the anchor provenance.
- **Orientation:** signed metrics (centroid shift, flank asymmetry) are defined only for
  **oriented** anchors (TSS, stranded motifs). For unoriented peaks, report unsigned or
  symmetrised versions and mark signed metrics as not applicable.
- **Strand-orient** oriented regions: reverse minus-strand rows and mirror the window around the
  anchor base, as deepTools/bwtool do. Preserve genomic coordinates and the mapping from bin to bp.
- **Grid:** primary ±2 kb with 20 bp bins (configurable); sensitivity at ±1/2/5 kb and 10/20/50 bp.
  No hidden smoothing is added, but note that bigWigs may already be smoothed (`bigwig_params`).
- **Fragment length and smoothing confound width.** Different fragment-size distributions
  between conditions, or read extension, produce apparent width changes. Compare `fragment_length`
  across conditions and report any difference. Where BAMs exist, regenerate coverage with
  identical parameters. For ATAC and CUT&RUN, prefer cut-site (5′-end, Tn5-shifted) coverage.
- **Fixed region set** for all conditions. If regions were discovered from the same data,
  document the selection bias and use union/consensus or matched controls.
- **Record:** library size, normalisation (CPM/RPGC/spike-in; ZEN is an optional normalisation
  comparison, not a shape test), coverage distribution, blacklist, signal sign, zero fraction and
  missingness (bins beyond chromosome ends are NaN, not 0).
- Signed or negative coverage is never passed into probability-distribution metrics without an
  explicit, reported transformation. Nothing is clipped silently.

## Module 1 — QC and amplitude (v1)
- **Per sample and region:** AUC (signal × bp), mean, maximum, centre and flank integrals,
  centre/flank ratio, missingness, zero fraction.
- **Counts from BAM (v2):** aggregate raw integer counts in non-overlapping or appropriately
  handled intervals and use DESeq2/edgeR with replicate/batch design.
- **Normalised bigWigs:** replicate-level comparisons of summaries. Never run DESeq2 on
  floating-point normalised signal.
- Show replicates separately. PCA/correlation are descriptive, not tests.

## Module 2 — shape / transport (v1)
**Background correction before normalisation.** If background is left in, a peak that loses height
over a constant background looks flatter, so W1 and width rise and an amplitude change is
misread as broadening. W1 weights distant mass, so flank background dominates it.
- Estimate the local background per profile from the outer flanks (default: median of the
  outermost 500 bp on each side) and subtract it.
- Compute shape on a **core window** (default ±1 kb), not the whole extraction window.
- Negative values after subtraction follow an explicit policy: `truncate` (set to 0 **and report
  the truncated mass fraction**) or `exclude`.

**Condition-independent inclusion.** A profile set enters the shape analysis only if the
**pooled** (all samples) background-corrected core mass passes a threshold, and every sample has
positive mass. Report the exclusions per condition. Per-condition thresholds are not allowed:
they bias towards the stronger condition.

**Metrics.** For each included profile, `p_i = y_i / Σ y_i` on physical bin coordinates:
- centroid (bp), SD width (bp), IQR width (bp);
- flank asymmetry (right − left mass)/total — oriented anchors only;
- between distributions: W1 (bp), W2 (bp), Jensen–Shannon divergence (bits).

**Required unit tests:**
- identical profiles → 0;
- pure translation by 200 bp → W1 ≈ 200 bp when fully contained;
- positive rescaling → normalised W1 = 0;
- asymmetry sign reverses under orientation reversal;
- W1 is correct with non-uniform bin widths.

**Levels of comparison.** Distinguish W1 between **group-mean distributions** from the
distribution of **region-level** distances. The latter compares condition means per region:
regions pair WT and KO, but replicates are not paired unless `pairing_id` says so. Report
amplitude and shape in separate panels. A naive version without background correction is kept as
a comparator only, to demonstrate the confound.

### v1 implementation findings (supersede the bullets above where they differ)
The pilot benchmarks showed that the background step needs four further safeguards. Without each
one, *amplitude-only change over a constant background* was still reported as a width change.
1. **Use the mean, not the median, of the flanks.** Count noise is right-skewed, so the median sits
   below the true background. The positive floor left after subtraction made taller peaks look
   narrower (false SD-change rate 37 % with the median vs 3 % with the mean).
2. **Do not truncate per region.** Clipping negative noise to zero recreates a positive floor.
   Shape is measured on each replicate's **aggregate** profile using **signed** moments. Centroid,
   variance, CDF/W1 and asymmetry are well defined for a positive-total signed measure, and
   zero-mean noise averages out. Only JSD needs non-negative input; it uses the clipped profile and
   is descriptive.
3. **Bound the per-region weights.** Each region is scaled by *pooled mass + median pooled mass*,
   which is condition-independent. Dividing by pooled mass alone is unstable at low SNR (near-zero
   denominators), and on real low-SNR TSS templates it produced CI coverage of 0.10–0.40.
4. **Benchmark truth must match the region set.** True effects are computed on the noise-free
   version of *each run's own* regions, with mean-one noise factors. A single shared truth
   understated coverage.

With these, both a parametric and a real-template benchmark (40 reps each) gave:
- false-call rate ≈ 0.03 (median) where the truth is ≈ 0;
- CI coverage ≥ 0.88 (median 0.97);
- amplitude-only shape false calls ≤ 0.10, while the naive metric was fooled 100 % of the time.

## Status (October 2026)
- **v1.1 done:** primary / secondary (Holm) / comparator roles, `log2_enrichment`, diagnostics
  (background ratio, amplitude heterogeneity / re-weighting, bidirectional shifts, replicate side of
  primary effects), per-sample bounded weighting (`shrunk_sample`), family-wise false-call rate in the
  benchmark.
- **v2 done (re-prioritised ahead of segmentation):** per-region empirical-Bayes moderated tests
  (`regional.py`, limma-identical), BH FDR, background-corrected relative amplitude with robust
  offsets, median-of-regions robust global shape, per-region FDR/power benchmark, `power` command.
- **Validation and comparison:** `benchmarks/profileseg_validation.py` (frozen method, criteria fixed
  in advance), `benchmarks/profileseg_comparison.py` (pointwise baselines, regions-as-units, edgeR-TMM).
- **Next (v3):** segmentation / interval-wise local inference (this module), BAM count path, robust
  empirical-Bayes hyperparameters (Phipson et al. 2016) for heavy-tailed per-region shape metrics.

## Module 3 — segmentation (v3; was v2)
- Contrast curve at aligned bins (raw difference and log-ratio with documented pseudocounts);
  penalised piecewise-constant segmentation (FPOP-style), with configurable penalty and minimum
  segment length. Record the library and cost function.
- **Where BAMs exist, use a DiffSegR-like count path.** DiffSegR segments a coverage log2FC and
  then fits DESeq2 to counts per segment *(verify)*. Inference comes from the count model, not
  the segmentation, and is still selective when the same data are reused.
- Labels: `Higher`, `Lower`, `Not_significant`, and `Equivalent` only with a pre-specified
  equivalence margin. Never call p > 0.05 "Same". Distinguish estimated boundaries from
  significant intervals; run penalty sensitivity.

## Module 4 — inference and multiplicity
**Two different targets; do not mix them.**
1. **Global** ("on average, does shape/amplitude change between conditions?"). The unit is the
   biological replicate. Per replicate, summarise each metric across regions, then contrast the
   conditions. No per-peak pseudoreplication.
2. **Per-region** ("which regions change?"). Test each region using replicates as observations,
   as csaw/DiffBind do, with multiplicity control across regions. This is legitimate and is not
   pseudoreplication. Count models apply when counts are available.

**Small-n reality.**
- Exact replicate-label permutation for 3 vs 3 has C(6,3) = 20 assignments. A symmetric
  (two-sided) statistic takes at most 10 distinct values, so **the smallest attainable two-sided
  p is 0.1, and 0.05 is unreachable.** For 2 vs 2 the floor is 1/3.
- Always report the attainable p-value floor. Never imply fine-grained p-values from many
  resamples of a tiny design.
- Credible inference at n = 2–3 needs **information borrowed across regions**: a hierarchical or
  mixed model (replicate as unit, region as random effect) or limma-style moderated statistics.
  If neither applies, report effect sizes with uncertainty and label the result exploratory.

**Uncertainty.** Bootstrap **hierarchically**: resample replicates within condition, then regions
(keeping regions paired across conditions). Bootstrapping regions alone understates the
uncertainty. With 2–3 replicates, bootstrap intervals are crude; say so.

**IWTomics** compares two populations of curves and treats **regions** as the units *(verify)*.
Its p-values therefore address region populations, which is in tension with replicate-level
inference. Apply it only after averaging replicates per region, and describe the target
accordingly.

**Selective inference.** Segmenting and then testing on the same data is selective. Compare
pre-specified intervals, sample splitting, and rerunning segmentation inside every null
resample. Correct across positions, segments, region sets, contrasts and metrics. BH alone does
not fix selection or dependence.

## Module 5 — benchmarks with known truth
**Semi-synthetic first.** Pure parametric simulation risks circularity: methods look good under
their own assumptions. So the benchmark uses:
- **Template-based simulation:** profile shapes taken from real data, with replicate noise
  (negative binomial, replicate-level amplitude dispersion, region heterogeneity, autocorrelation).
- **Replicate-split null:** split ≥ 4 real replicates of one condition into two pseudo-groups, to
  calibrate false positives on real variability.
- **Injection:** apply known transformations (shift, broadening, amplitude, flank) to real
  profiles.

**Scenarios** (at least):
1. Null.
2. Uniform amplitude gain/loss **of the peak over a constant background** (key confound test).
3. Global multiplicative scaling including background (erased by normalisation).
4. Pure translation.
5. Broadening/narrowing at conserved area.
6. Left-flank-only and right-flank-only changes.
7. Bimodal/shoulder change.
8. Sparse subset of affected regions.
9. Low coverage, zero inflation, misalignment, library-size change.
10. **Fragment-length difference** between conditions (apparent width).
11. **Anchors chosen from one condition** (apparent shift).

**Methods:** AUC contrasts; background-corrected shape metrics; **naive uncorrected** shape
metrics; (v2) FPOP-only segmentation, the DiffSegR-like count path and IWTomics. DiffSegR is run
only on supported inputs.

**Report** null false-positive rate, power, localisation error (bp), direction accuracy,
effect-size bias, CI coverage and runtime/memory, with fixed seeds, repetitions and denominators.
Never invent results; pilot honestly and queue full runs.

## Module 6 — publication figures
Vector PDF/SVG + PNG:
- A: per-replicate profiles with group means and uncertainty;
- B: raw vs background-corrected, area-normalised profiles;
- C: difference curve with uncertainty band and zero line (no significance claimed in v1);
- D (v2): significant intervals distinct from FPOP boundaries;
- E: effect sizes (ΔAUC, Δcentroid, Δwidth, W1, asymmetry) with CIs;
- F: benchmark error/power/bias.

Reuse GenomeCanvas for extraction, which is validated against bwtool and handles minus strands
and chromosome edges, and for styling. Colours, fonts, sizes and line widths are configurable.

## Engineering and reproducibility
- Python package `genomecanvas.profileseg` with CLI `gcanvas-profileseg` and a YAML config; pinned
  environment. Prefer validated libraries for complex inference.
- **Unit tests:** indexing, strand reversal, bins → bp, all-zero/negative input, empty regions,
  ties, NaN, non-uniform bins, overlapping regions, permutation enumeration.
- Synthetic fixtures and a small end-to-end example runnable on a laptop.
- Log commands, parameters, versions, seeds, exclusions and runtime. Produce human-readable
  reports plus TSV outputs.
- Stream/chunk extraction; submit cluster (LSF) jobs only after confirming policy and resources.

## Execution gates
- **Gate 1:** discovery, manifest, package feasibility, design report. Stop and ask if replicates,
  anchors, genome build or contrast cannot be established.
- **Gate 2 (v1):** prototype on synthetic data with tests, plus template-based pilot benchmark.
- **Gate 3:** real-data exploratory run and QC.
- **Gate 4 (v2):** segmentation, local inference, calibration; report failures too.

## Acceptance criteria (v1)
- A single documented command runs the synthetic example end-to-end, including bigWig extraction.
- Amplitude-only, shift-only, width-only and flank-only cases are each described correctly by the
  effect sizes. In particular, **amplitude-only over a constant background is not reported as a
  shape change** after background correction, while the naive metrics demonstrate the confound.
- The attainable permutation p-value floor is reported for the design actually used.
- No significance claims from segmentation alone. Replicates are never silently treated as
  independent regions.
- Reports separate **observed**, **inferred** and **hypothesised** conclusions, and include
  limitations.

## Final agent response
Summarise:
1. discovered data;
2. exact scripts and paths;
3. implemented vs pending methods;
4. tests and benchmark results with denominators;
5. statistical limitations;
6. reproduction commands;
7. decisions needed from the user.

Do not report unexecuted steps as completed.
