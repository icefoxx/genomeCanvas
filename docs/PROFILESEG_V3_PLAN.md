# ProfileSeg v3 — plan

Baseline: tag `profileseg-v2.4` (commit `13206cc`). Development branch: `profileseg-v3`.

## Goals

1. **Two first-class input paths, adapted and evaluated separately.**
   - **bigWig** (coverage tracks as delivered by most pipelines and GEO): unknown scaling, smoothed
     coverage, no fragment information.
   - **BAM** (aligned reads): integer fragment counts, fragment lengths, strand, duplicates and
     mapping quality are available, so the method can use them.
2. **Close the open v2.4 limitations** that the input type can address:
   per-region absolute amplitude (miscalibrated in some designs), no depth diagnostic (sparse data
   gave no power without warning), aggregate width metrics under amplitude heterogeneity.
3. Keep the v2.4 per-region shape tests unchanged unless an evaluation shows a gain.

## What each path provides

| capability | bigWig | BAM |
|---|---|---|
| signal per bin | coverage as stored (scale unknown) | fragment midpoints, extended coverage, or read ends (strand-aware), from integer counts |
| depth diagnostic | *estimated* fragments per region from the track's count quantum (smallest positive increment) when it exists; otherwise relative only | exact fragments per region and sample |
| per-region amplitude | moderated t on log signal (v2.4) | counts with precision weights (voom-style log-CPM + weighted moderated t) |
| fragment-length QC | not available | per-sample fragment-length distribution (paired-end), flags a deviant replicate |
| fragment-size selection | not available | e.g. nucleosomal vs sub-nucleosomal fragments (MNase, CUT&RUN, ATAC) |
| filters | as produced upstream | MAPQ, proper pairs, duplicates, size range |

## Evaluation (bigWig vs BAM, separately)

1. **Simulation with both outputs from the same counts**: BAM (reads placed from simulated
   fragment counts) and a coverage track from the same fragments; run both paths; compare
   calibration (nulls, amplitude-only changes) and power (shifts, broadening, site amplitude).
2. **Real data, same reads**: yeast MNase and TF CUT&RUN BAMs, and coverage tracks built from the same
   BAMs, through both paths; real nulls and read-level insertion.
3. **Depth diagnostic**: does it predict where per-region shape tests have power (yeast ~400 fragments
   per region vs CUT&RUN ~10)?
4. Frozen validation with pre-specified criteria on fresh seeds for each path before release.

## Order of work

1. BAM reader (`bamio.py`: midpoint / coverage / end modes; MAPQ, proper-pair, duplicate and size
   filters; paired- and single-end) and manifest support (`.bam` signal files).
2. Depth diagnostic for both paths.
3. Count-based per-region amplitude for BAM (voom-style weights); evaluate against v2.4 on the designs
   where absolute amplitude failed.
4. Fragment-length QC and size selection.
5. bigWig vs BAM evaluation (simulation, then real data).
6. Aggregate robust estimators (median-of-regions shape metrics as the reported aggregate).

## Progress (2026-10-26)

| item | status | result |
|---|---|---|
| BAM reader (`bamio.py`) | done | midpoint / coverage / end5 / end3; MAPQ, proper pairs, duplicates, size range; paired-end via TLEN; single-end extension with strand cross-correlation estimate; library size, per-region fragments, fragment-length summary. Cross-checked against samtools-based counts on yeast BAMs (totals equal; 1 bp coordinate convention difference) |
| manifest support | done | `.bam`/`.cram` signal files select the BAM path; `bam.*` options in the config |
| depth diagnostic | done | BAM: exact fragments per region; bigWig: count quantum + single-fragment footprint read from the track (quantum found in 4/4 GEO tracks; footprints 2 bp point, 60 bp read coverage, 157 / 190 bp fragment coverage); yeast estimate ~1,000–1,450 vs ~1,600 expected (within ~1.5×). Warning below a median of 50 fragments |
| count-based amplitude (voom weights) | **evaluated, rejected** | calibrated on Poisson-only counts but not on ChIP-like regional counts (false share 1.4–6.7 % of exact-null regions): variance there is region-specific biology, not count-driven |
| exact-null per-region evaluation | done | amplitude tests: ≤ 1.3 % false calls among exact-null regions, but complete-null FWER-type R2 still 0.20 (3 vs 2) and 0.275 (2-cell templates) for absolute `log2_signal`: a real, mild liberalness remains |
| bigWig vs BAM (simulation, same fragments) | done | both calibrated (0 calls under null / amplitude-only); power similar, coverage (bigWig) slightly higher than BAM midpoints (smoother); Jaccard agreement 0.80–0.99 |

Next: BAM coverage vs midpoint mode on real data; bigWig vs BAM on the same real reads (yeast, CUT&RUN);
calibrate the depth threshold by simulation; size selection; aggregate robust estimators; a fix for the
residual liberalness of per-region absolute amplitude (e.g. robust eBayes, Phipson 2016).
