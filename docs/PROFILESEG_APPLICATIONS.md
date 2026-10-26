# ProfileSeg — where to use it, and how it compares with other tools

## 1. What ProfileSeg answers

Given **predefined anchors** (TSS, ATG, peak summits, motif sites, nucleosome dyads, cut sites)
and **coverage tracks with biological replicates** in two conditions, ProfileSeg asks, separately:

1. **Amplitude**: is there more or less signal (aggregate and per region; absolute, relative to
   the bulk of regions, or relative to the local background)?
2. **Shape**: did the signal move, spread or change form (aggregate and per region), independently
   of amplitude?

Replicates are always the unit of inference. Per-region shape calls come from `shape_kernel`
(omnibus) and `shift` (directional, with bp effect size); both use a null pooled over regions and
replicate-label assignments, so they work with 2 vs 2.

## 2. Application scenarios

| assay / question | anchors | what ProfileSeg gives | status of evidence |
|---|---|---|---|
| **MNase-seq / ATAC nucleosome organisation** (remodeler mutants, drugs): +1 nucleosome shift, loss of phasing, fuzziness | ATG/TSS (oriented) or nucleosome dyads | which genes change array shape; direction and size of shifts; amplitude (occupancy) kept apart | **real data**: yeast isw1Δ chd1Δ (GSE280671), see `PROFILESEG_COMPARISON.md` §6 |
| **Pol II ChIP / PRO-seq / PRO-cap**: pause-site position, pause release, initiation focus | TSS (oriented), strand-matched tracks | pause-peak shift (bp), broadening, promoter-proximal vs gene-body signal | real-data check with NELF degradation (GSE233930), §6 |
| **TF ChIP / CUT&RUN / ChIP-exo**: binding loss vs redistribution; mutant TFs that change where they sit relative to nucleosomes | peak summits, motif sites, nucleosome dyads | amplitude change per site, and whether binding moves (e.g. toward linker / NFR) | simulation + semi-real |
| **DNA damage / cut-site profiles** (restriction-enzyme or Cas9 sites; resection, repair-factor spreading) | cut sites | signal level vs spread around the break, per site | simulation |
| **Histone marks around TSS / enhancers** (H3K4me3 domain width, H2A.Z, acetylation spreading) | TSS, enhancers | domain spreading or shift vs level change | simulation; broad domains (> 5 kb) are out of scope |
| **Comparing tagged constructs or antibodies** (Flag vs HA, mutant vs WT protein) | binding sites | whether apparent differences are amplitude, background or shape | simulation |

Design guidance (from validation, §4–5 of `PROFILESEG_COMPARISON.md`):

| replicates per group | what can be called per region | aggregate (typical-region) statements |
|---|---|---|
| 1 | nothing (descriptive only; region-level p-values would be pseudo-replication) | descriptive |
| 2 | shifts ≥ ~100 bp (power 0.5–0.9); not 50 bp; broadening rarely | amplitude ≥ 1.2×, shifts ≥ 25 bp (power ~0.75) |
| 3 | shifts of 100 bp (0.8–0.97), 50 bp in 15–45 % of regions; broadening 0.3–0.6 | most effects of interest |
| ≥ 4 | gains roughly √n in z; best investment for small per-site effects | — |

Not suitable:
- genome-wide discovery of differential regions without anchors (use csaw / THOR / DiffBind);
- very broad domains (H3K27me3, H3K9me3 blocks): the window model assumes a local profile;
- translational efficiency or other ratio-of-assays questions (Ribo-seq TE: Xtail, deltaTE);
- **absolute** global amplitude claims without spike-ins (ProfileSeg flags but cannot resolve them).

## 3. Comparison with tools of the same kind

Claims about other tools come from their publications and, where stated, from our own runs
(`PROFILESEG_COMPARISON.md`). Tools marked "not run here" were not benchmarked by us.

| tool | question | input | replicates as unit | amplitude vs shape | per-region shape | effect sizes | works at 2 vs 2 | our evidence |
|---|---|---|---|---|---|---|---|---|
| **ProfileSeg** | amplitude and shape, aggregate + per region | bigWig / computeMatrix | yes | separated | yes (`shape_kernel`, `shift`) | log2FC, bp | yes | validated on simulation, semi-real noise and real yeast MNase nulls / positives |
| **MMDiff2** (Schweikert 2013) | read-distribution difference per peak | BAM + peaks | null from replicate comparisons | not separated (responds to amplitude) | yes (omnibus) | no | partly | run here: 34 % of regions called under the null with its own p-values (Fisher combination of dependent pairwise p-values); calibrated, lower power than `shape_kernel` |
| **IWTomics** (Cremona 2018) | where along the window two **region sets** differ | aligned feature matrices | no (regions are curves) | not separated | no (localises within the window) | no | no (3 vs 3 gives p ≥ 0.1 with replicates as curves) | run here: 100 % calls under amplitude-only changes with regions as units |
| **DANPOS2/3** (Chen 2013) | per-nucleosome position / fuzziness / occupancy change | MNase BAM/wig | Poisson-type tests on (pooled) reads, as far as we know | separated by decomposition | per nucleosome | yes | not documented | not run here |
| **Nucleosome Dynamics** (Buitrago 2019) | nucleosome shifts, evictions, fuzziness changes | MNase reads | not as replicate-level inference, as far as we know | separated | per nucleosome | yes | — | not run here |
| **nucMACC** (preprint 2022) | nucleosome accessibility / stability from MNase titrations | MNase titrations | — | accessibility, not shape tests | — | — | — | not run here |
| **DiffBind** (Ross-Innes 2012) / **csaw** (Lun & Smyth 2016) | which peaks / windows change in count | BAM | yes (edgeR / DESeq2) | amplitude only | no | log2FC | yes | edgeR-TMM run here: same power as ProfileSeg relative amplitude; slightly anti-conservative with background-containing windows |
| **THOR** (Allhoff 2016) | differential peaks via HMM with replicates | BAM | yes | amplitude only | no | — | yes | not run here |
| **MAnorm2** (Tu 2021) | differential signal with hierarchical normalisation | counts | yes (eBayes) | amplitude only | no | log2FC | yes | not run here |
| **PePr** (Zhang 2014) | differential peaks with replicates | BAM | yes | amplitude; optional shift-size filter | no | — | yes | not run here |
| **deepTools** plotProfile / plotHeatmap | visualisation | bigWig | no inference | — | — | — | — | — |

### Where ProfileSeg is the better choice

- You have **coverage tracks** (not necessarily BAMs) and replicates, and the question is whether
  binding **moved or spread** rather than only went up or down.
- You need to **separate amplitude from shape** (e.g. a mutant that binds less *and* elsewhere).
- You have **two replicates per group** and need calibrated per-region shape calls.
- You want **effect sizes** (bp, log2) and diagnostics (background change, re-weighting, replicate
  artefacts) rather than a single distance or a significance map.

### Where another tool is better

- **Count differential binding with BAMs**: DiffBind / csaw (mature count models, spike-ins).
- **Per-nucleosome calls genome-wide** (without predefined anchors): DANPOS / Nucleosome Dynamics;
  ProfileSeg can then test the called dyads with replicates.
- **Localising where along a window two different region sets differ**: IWTomics.
- **Any distributional difference at a peak, including complex shapes, with BAMs**: MMDiff2 as a
  ranking tool (its p-values need external calibration and an amplitude check).

## 4. Known limitations

- Anchors must be given; windows are local (± a few kb).
- Shape tests need predominantly positive background-corrected signal (on MNase data about half of
  the yeast ORFs qualify with the default window; widen the core or choose dyad anchors).
- Small per-site effects are information-limited (2 vs 2: 50 bp shifts are not detectable).
- The template-fit shift estimate is attenuated by up to ~27 % at 2 vs 2 on broad profiles.
- Replicate heterogeneity that aligns with the contrast (e.g. a severity gradient among mutant
  replicates) is detected as a difference — correctly, but it means replicates must be
  representative. Sample-wide artefacts in a single replicate are absorbed by the pooled null and
  flagged by the replicate-influence diagnostics.
