# ProfileSeg — manuscript plan (working draft)

## Positioning

**Core claim (narrow and defensible).** With as few as two biological replicates per condition,
ProfileSeg gives calibrated per-region tests of whether the *shape* of a coverage profile changed
(omnibus `shape_kernel`; directional `shift` with an effect size in bp), separately from changes in
*amplitude*, directly on coverage tracks.

**Secondary claims.**
1. Existing approaches are miscalibrated in replicated designs for identifiable reasons:
   MMDiff2 combines dependent pairwise p-values with Fisher's method (34 % of regions called under
   the null in our simulation); region-as-unit methods (IWTomics in its intended use, Wilcoxon over
   peaks) are pseudo-replicated and respond to amplitude changes.
2. A label-invariant stabilising constant is necessary for pooled-permutation shape tests: an s0
   computed per label assignment makes the test anti-conservative when many regions change.
3. Diagnostics for confounds that otherwise masquerade as shape changes (background change,
   amplitude re-weighting, a deviant replicate).

**Not claimed.** The aggregate (typical-region) analysis is presented as descriptive support with
explicit conditions of use (`docs/PROFILESEG.md`, "When to use the aggregate centroid"); its earlier
validation failures are reported, not hidden. No claim of genome-wide discovery without anchors.

## Evidence map

| claim | evidence in hand | still needed |
|---|---|---|
| calibration with 2–3 replicates | three frozen validations on fresh seeds with pre-specified criteria (V1–V3 pass in all cells); semi-real data with real replicate noise; real nulls: 18 yeast 2 vs 2 splits (`shift`: 0 calls in all; `shape_kernel`: 0 in 15, calls only where the split coincides with batch-like replicate pairs) and four human PRO-seq 3 vs 3 line-balanced splits (0 calls); insertion on real replicates (FDP 0.03–0.05 yeast, < 0.02 PRO-seq) — `PROFILESEG_COMPARISON.md` §8 | a 3 vs 3 real null from a single cell line; batch-aware label assignments |
| power and effect sizes | simulation, semi-real, insertion on real replicates; PRO-seq NELF: 98 % downstream, +17.5 bp vs ~18 bp independent meta-profile | — |
| separation of amplitude and shape | amplitude-only scenarios in every benchmark, including insertion on real data | — |
| competitor miscalibration | simulation (MMDiff2, IWTomics, edgeR); MMDiff2 source analysis; **real BAMs** (yeast MNase realigned from ENA reads): MMDiff2 called 7–17 of 134 genes in real nulls vs 0 for ProfileSeg, and 27 vs 65 in the positive control (§9); conventional practice on real data: Wilcoxon over regions significant in 78–100 % of real null splits, pausing index FWER 0.75 (§8) | a second assay with BAMs |
| practicality | bigWig input; runtime/memory vs MMDiff2 (simulated BAMs, 300–1000 regions: 3–8 s vs 7–15 min, 0.3–0.4 vs 6.6–7.8 GB); exact kernel feature map (19× faster default pipeline at 20,000 regions, 4 vs 4) — §10 | packaging as a stand-alone tool (planned separately) |

## Figures (draft)

1. Problem: amplitude vs shape confounds; why regions/reads are not replicates (schematic + IWTomics /
   MMDiff2 null behaviour).
2. Method: per-region statistics, pooled permutation null over regions and label assignments,
   label-invariant s0 (with the per-assignment s0 failure as a panel).
3. Calibration: frozen validation grid (sources × designs), real nulls, insertion on real data.
4. Power: shift / broadening vs replicates; small-shift information limit (realised z).
5. Real data: yeast remodeler mutants (phasing check), PRO-seq NELF (direction and size).
6. Comparison and resources: MMDiff2, IWTomics, edgeR, conventional practice (Welch on centroids,
   Wilcoxon over regions, pausing index); runtime and memory.

## Threats to validity (to state in the paper)

- Simulator and tuning seeds share a generative family; mitigated by fresh-seed validations,
  semi-real noise and insertion on real replicates, but not eliminated.
- Development history: the selection rule was amended (FDP constraint) and two defects were fixed
  (per-assignment s0; dependent semi-real noise) before the final validations; all are documented.
- The `shift` effect size is attenuated (up to ~27 %) at 2 vs 2 on broad profiles (criterion V4).
- Replicate heterogeneity aligned with a contrast is detected as a difference (yeast double mutant
  split); interpretation relies on representative replicates.
- Small per-site effects are information-limited (50 bp at 2 vs 2).

## Target venues

Mid-tier methods journals (NAR Genomics and Bioinformatics, BMC Bioinformatics, Bioinformatics)
once the in-progress items are complete and the tool is packaged on its own.
