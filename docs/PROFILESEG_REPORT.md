# ProfileSeg evaluation report

*Name: "Seg" currently means region-level screening of differential profiles, not segmentation of curves into changed intervals.*

*Status as of 2026-10-26, end of day (code: ProfileSeg v2.4, uncommitted working tree after commit `a6fa276`).
All numbers come from the scripts in `benchmarks/` and are detailed in `docs/PROFILESEG_COMPARISON.md`.*

## 1. Summary

ProfileSeg compares coverage profiles around predefined anchors between two conditions with
biological replicates. It reports **amplitude** (how much signal) and **shape** (where and how
spread the signal is) separately, for the typical region and for each region, and treats replicates
— not regions or reads — as the unit of inference.

**Main result.** The per-region shape tests (`shape_kernel`, omnibus; `shift`, directional with an
effect size in bp) are calibrated with as few as **two replicates per condition** on simulated,
semi-real and real data, and they detect known changes in two real assays. On the same real BAMs,
the closest existing tool (MMDiff2) made false calls in real nulls where ProfileSeg made none, and
detected fewer genes in a known positive control.

**Main limitations.** Small per-site effects are information-limited (50 bp shifts at 2 vs 2);
processing batches that coincide with the contrast are detected as differences; the `shift` effect
size is attenuated by up to ~27 % at 2 vs 2 on broad profiles; the aggregate (typical-region)
analysis is reliable only under stated conditions. Evidence for some application areas (TF
redistribution, cut sites, histone spreading) is simulation-only.

| aspect | assessment |
|---|---|
| statistical design | replicate-level inference; shape metrics use background-corrected, mass-normalised profiles, which separates shape from amplitude exactly for noise-free multiplicative changes and approximately with noise and imperfect background estimates (checked by amplitude-only scenarios, not guaranteed); per-region p-values come from a null pooled over regions and label assignments, so they depend on the null being comparable across regions (§4, §7) |
| calibration evidence | strong for per-region tests: frozen validations (final round passes every criterion), real nulls in 3 assays (yeast MNase 2v2, PRO-seq 3v3, single-line CUT&RUN 3v3), insertion on real replicates |
| power | per region (20 % of regions changed): 100 bp shifts 0.52–0.97; 50 bp shifts 0.16–0.47 at 3 vs 3 and ≤ 0.16 at 2 vs 2; partial broadening 0.01–0.65 depending on shape and design; coherent ~20 bp shifts detected on real PRO-seq at 2 vs 2 |
| robustness | single-replicate artefacts absorbed and flagged; batch structure aligned with the contrast is detected (correctly, but it must be designed out) |
| practicality | bigWig / computeMatrix input; seconds for 5,000 regions; ~100× faster and ~20× less memory than MMDiff2 |
| maturity | beta; embedded in GenomeCanvas; single developer; not yet packaged separately |

## 2. Method in brief

| component | what it does |
|---|---|
| background correction | per sample and region, mean of the outer flanks |
| amplitude (per region) | `log2_signal`, `log2_signal_relative` (relative to the bulk of regions, half-sample-mode offsets), `log2_enrichment`; limma-identical moderated t, BH |
| `shape_kernel` | Gaussian-kernel (MMD-type) distance between condition-mean mass-normalised profiles, divided by within-condition spread + s0; null pooled over regions × non-observed label assignments |
| `shift` | projection of each sample's profile on the derivative of the pooled template (locally most powerful for small translations), same pooled null; effect size in bp from a weighted Gauss–Newton template fit |
| label-invariant s0 | s0 from the total spread, identical under every assignment (a per-assignment s0 made the test anti-conservative; found and fixed during development) |
| diagnostics | background ratio, amplitude re-weighting, replicate-wide offsets relative to the group, leave-one-replicate-out retention (≥ 3 per group) |
| options | condition-blind anchor alignment (FunChIP-style shift-only k-mean alignment), descriptive shape classes |

## 3. Evidence

### 3.1 Simulation and semi-real data (frozen validations)

Criteria fixed before running: V1 null FWER ≤ 0.10; V2 false-call share ≤ 0.005 in no-shape
scenarios; V3 FDP ≤ 0.10 in 20 % mixtures; V4 shift estimate within ±20 %.

| round | seeds | what changed | result |
|---|---|---|---|
| 2 | 700 000+ | shift score test | V1–V3 pass; V4 fails (real templates 2v2) |
| 4 | 2 500 000+ | batch handling, refined shift estimator | all pass |
| 5 | 3 000 000+ | kernel: whole-window signed normalisation | calibrated, but kernel power collapsed on broad profiles (rejected) |
| **6 (final)** | 3 500 000+ | kernel: signed profile / support mass | **all criteria pass in all 6 cells**; zero false calls under amplitude-only changes |

Final per-region power (round 6), kernel / shift at 100 bp: 0.52–0.91 / 0.56–0.97; shift at 50 bp
0.01–0.47; `shift` estimate 84.5–99.9 bp for a 100 bp truth.

Global validation (round 4, revised criteria): per-region shape tests pass in all six designs;
aggregate width metrics fail A2/A3 under amplitude changes in some designs (re-weighting confound),
and per-region *absolute* amplitude fails R2 in three designs (§4).

### 3.2 Real data

| dataset | design | result |
|---|---|---|
| yeast MNase, WT replicate splits (GSE280671) | 2 vs 2, 3 splits | 0 calls (both tests) |
| yeast MNase, 18 replicate splits over 6 genotypes | 2 vs 2 | `shift`: 0 calls in all; `shape_kernel`: 0 in 15; calls (5–18 % of regions) only in "replicates 1+2 vs 3+4" splits of three genotypes whose replicate pairs differ systematically (cross-pair / within-pair distance 1.5–2.7) |
| human PRO-seq DMSO, line-balanced (GSE233930) | 3 vs 3, 4 splits | 0 calls (both tests) |
| insertion on real replicates (20 % of regions) | yeast 2 vs 2, PRO-seq 3 vs 3 | FDP 0.11–0.13 with all yeast splits and 0.03–0.05 after excluding three batch-aligned splits (exclusion rule defined after seeing the results, i.e. post hoc), < 0.02 (PRO-seq); `shift` 20 bp power 0.97–1.00; `shape_kernel` 100 bp 0.99–1.00 |
| yeast isw1Δ chd1Δ vs WT | 4 vs 4 | 91 % of genes change shape; ordering of mutants agrees with an independent phasing measure (WT 0.97, singles 0.89–0.94, double 0.40–0.67) |
| human PRO-seq, NELFB degradation | 2 vs 2 | 6,350 genes shifted, 98 % downstream, +17.5 bp (independent meta-profile ~+18 bp) |
| yeast MNase, same real BAMs, vs MMDiff2 | 2 vs 2, 134–141 genes | nulls: ProfileSeg 0 and 0, MMDiff2 33 and 20 (of 300); positive: ProfileSeg 90 of 141 |
| TF CUT&RUN, single cell line (GSE291134) | 3 vs 3, 2 lines × 10 splits | per-region and aggregate: 0 calls (one split 0.14 %); but ~4–20 fragments per site: no power, MMDiff2 crashed — the second real-BAM comparison is inconclusive |

### 3.3 Conventional practice on real data

| practice | behaviour |
|---|---|
| Wilcoxon over regions (regions as units) | significant in 78 % (yeast) and 100 % (PRO-seq) of real null splits |
| Welch t on window centroids | no power on MNase arrays; low on PRO-seq |
| moderated t on window centroids (ProfileSeg v2) | FWER 0.25 on PRO-seq nulls |
| pausing index (moderated t) | FWER 0.75 on PRO-seq nulls; ≤ 0.27 power for shifts |

### 3.4 Runtime and memory

Default pipeline (current code, single thread): 2.8–7.9 s for 5,000 regions (2v2–4v4), 11–36 s for
20,000. Peak memory 0.8–1.5 GB (5,000 regions) and 2.7–3.6 GB (20,000). An exact kernel feature map
made the default pipeline 19× faster at 20,000 regions (4 vs 4). MMDiff2 on simulated BAMs: 418 s /
7.8 GB (300 regions), 912 s / 6.6 GB (1,000 regions).

## 4. Overall assessment

**Strengths**

1. Addresses a real gap: replicate-level, per-region shape inference that is separated from
   amplitude, on coverage tracks, at 2–3 replicates.
2. Calibration is supported by pre-specified criteria on fresh seeds, real nulls in two assays and
   insertion on real replicates — more than most published tools in this area.
3. Explains why existing approaches fail (Fisher combination of dependent p-values in MMDiff2;
   pseudo-replication in region-as-unit tests; amplitude–position mixing in pausing indices) and
   shows it on real data.
4. Interpretable outputs (bp, log2), diagnostics, low resource use.

**Weaknesses and threats to validity**

1. Simulation-based tuning shares a generative family with the simulation validations (mitigated,
   not eliminated, by semi-real and real-data checks).
2. Development history includes an amended selection rule and two fixed defects (per-assignment s0;
   dependent semi-real noise) before the final validations.
3. Batch structure aligned with the contrast produces calls; only balanced designs (or batch-aware
   assignments, in progress) avoid this.
4. Information limits for small per-site effects at 2 replicates.
5. Effect-size attenuation of `shift` at 2 vs 2 on broad profiles.
6. Aggregate analysis: earlier validation failures, low power at 2 vs 2, unstable on periodic
   profiles; positioned as descriptive with explicit conditions of use.
7. Real-BAM comparison so far in one assay and at subsampled depth; MMDiff2 had to run in
   single-end mode (its paired-end mode fails in 1.30.0).

**Readiness.** Suitable for internal use and for reporting per-region shape changes with the
caveats above. For publication in a mid-tier methods journal: about 85 % ready for a paper centred on
the per-region shape tests; the open items are listed in §7.

## 5. Application scenarios

| scenario | anchors | question answered | evidence |
|---|---|---|---|
| chromatin remodelling (MNase-seq, ATAC) | ATG/TSS (oriented), dyads | which genes' nucleosome arrays shift or lose phasing; occupancy separately | real data (yeast) |
| transcription pausing (PRO-seq, PRO-cap, Pol II ChIP) | TSS, strand-matched | redistribution of promoter-proximal signal; a translation-model shift estimate (bp); release vs amount | real data (human) |
| TF mutants, cofactor loss (ChIP, CUT&RUN) | summits, motifs, dyads | binding lost vs moved (e.g. toward NFR) | simulation, semi-real |
| DNA damage / cut sites | cut sites | signal level vs spread around breaks | simulation |
| histone marks at TSS/enhancers | TSS, enhancers | domain spreading vs level change (local domains only) | simulation |
| construct / antibody comparisons | binding sites | amplitude, background or shape | simulation |

Design guidance: 1 replicate — descriptive only; 2 — per-region shifts ≥ ~20 bp when coherent,
≥ 100 bp in general; 3 — most effects of interest; ≥ 4 — small per-site effects. Balance processing
batches across conditions. Use spike-ins for absolute amplitude claims.

Not suited to: anchor-free discovery; very broad domains; ratio-of-assays questions (Ribo-seq TE);
absolute global changes without spike-ins.

## 6. Comparison with existing tools

| tool | question | replicates as unit | amplitude vs shape | per-region shape | effect sizes | 2 vs 2 | evidence here |
|---|---|---|---|---|---|---|---|
| **ProfileSeg** | amplitude + shape, aggregate and per region | yes | separated | yes | bp, log2 | yes | simulation, semi-real, real nulls, insertion, real positives, real BAMs |
| MMDiff2 (1.30.0, as run here) | read-distribution difference per peak | null from replicate pairs | not separated | omnibus | no | partly | in this version and our configuration: anti-conservative p-values (simulation 34 %; real yeast nulls 7–33 of 300 genes), paired-end mode fails, slow; this is evidence about the implementation and setting, not a refutation of the kernel idea, which ProfileSeg also uses |
| IWTomics | where along a window two **sets of curves** differ (localisation, all scales) | no (curves are regions in its intended use) | depends on input normalisation | localises within the window | no | no | answers a different question; on raw curves, rejecting under amplitude-only changes is correct behaviour, not a shape false positive; a same-target comparison (normalised inputs, replicate-level design) is still to do |
| DANPOS / Nucleosome Dynamics | per-nucleosome changes genome-wide | not replicate-level inference, as far as we know | separated | per nucleosome | yes | — | not run |
| DiffBind / csaw / edgeR | count differences (csaw: per window, aggregated to regions) | yes | DiffBind: amplitude; csaw's window tests can detect local redistribution but give no translation estimate | csaw: indirectly (windows) | log2FC | yes | edgeR as powerful as ProfileSeg's relative amplitude; csaw not yet compared on the same target |
| conventional practice (Wilcoxon over regions, window centroids, pausing index) | various | no / partly | mixed | no | — | — | miscalibrated or underpowered on real data (§3.3) |

## 7. Gaps: closed and open

Closed in this round: batch-aware label assignments (confounded designs refused; optional
within-batch permutation); `shift` attenuation (V4 passes everywhere); a single-cell-line 3 vs 3 real
null (CUT&RUN, 20 splits); kernel bias under amplitude changes (found by validation, fixed); a
miscalibrated exploratory test (`log2_scale`) removed; 19× faster pipeline.

Open:
1. Aggregate width metrics under amplitude heterogeneity (A2/A3) and per-region absolute amplitude
   in some designs.
2. A second real-BAM comparison with adequate depth (the CUT&RUN data were too sparse).
3. Batch information for public data (inferred here, not recorded in GEO).
4. Stand-alone packaging.
