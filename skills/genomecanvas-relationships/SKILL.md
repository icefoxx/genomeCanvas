---
name: genomecanvas-relationships
description: GenomeCanvas sample-relationship figures - correlation heatmap, PCA scatter (with interactive HTML) and PairGrid scatter matrix over BED regions, genome-wide bins or sliding windows (like deepTools multiBigwigSummary + plotCorrelation/plotPCA), plus Sample QC that flags mislabelled, swapped and low-quality samples in large collections from a sample sheet and a bins x samples matrix.
---

# Sample relationships

Run with the entry skill's runner ([genomecanvas-plot](../genomecanvas-plot/SKILL.md)).

| Plot slug | Draws | Extra outputs |
|---|---|---|
| `relationships/correlation` | clustered correlation heatmap (Ward) | `summary.tsv`, `figure.tsv` |
| `relationships/pca` | PC1 vs PC2, labels repelled, group markers | `pca_scores.tsv`, `pca.json`, `pca_interactive.html`, `summary.tsv` |
| `relationships/pairgrid` | pairwise scatter matrix with histograms | `summary.tsv` |
| `relationships/sample-qc` | QC overview (mark PCA, mark check, status per mark) | see [Sample QC](#sample-qc-large-collections) |

The first three build a sample × feature table (`summary.tsv`), then plot it.

## Step 1: choose what to summarise over (`rel_mode`)

| `rel_mode` | Features | Needs a BED? |
|---|---|---|
| `regions` (default) | each interval of the ticked region files | **yes**, ≥1 region file |
| `bins` | non-overlapping genome bins of `rel_bin` bp | **no**: region files are ignored even if listed |
| `windows` | sliding windows of `rel_bin` bp every `rel_step` bp | **no** |

In `bins` / `windows` mode the only region file that is used is the optional **blacklist**. Bins
overlapping it are removed. It does not need to be listed in `files`.

## Step 2: settings

| Key | Default | Applies to | Notes |
|---|---|---|---|
| `rel_bin` | 10000 | bins, windows | 10 bp – 100 Mb |
| `rel_step` | `rel_bin / 5` | windows | `rel_bin` must be a multiple of the step |
| `rel_scope` | `genome` | bins, windows | `chroms`, or `locus` (uses `chrom/start/end`) |
| `rel_chroms` | | scope `chroms` | `"chr1, chr2, chrX"`; `chr` prefix optional |
| `rel_skip` | `chrM, MT, *_random, chrUn*, *_alt, *_fix, chrEBV, *hap*, *decoy*` | bins, windows | case-insensitive wildcards |
| `rel_blacklist` | none | bins, windows | local: a BED path; server: file id / name / label |
| `rel_stat` | `sum` for regions, `mean` for bins/windows | all | `mean` (covered bases), `mean0` (missing = 0), `sum` (signal × bp), `max` |
| `rel_skip_zeros` | true | all | drop features that are zero in every sample |
| `rel_max_rows` | all | all | reproducible random subset (seed 7201) for speed |

Only chromosomes present in **every** signal are used, each trimmed to the shortest length among the
files. Limits: 3,000,000 rows after filtering, and about 15 M windows × samples before filtering.
Use larger bins, a larger step, `rel_chroms` or `rel_max_rows` to stay under them.

Speed: for genome-wide bins/windows and large BED sets, each sample's BigWig is read in its own
process (up to 4 at once), with identical results to a serial read. Six mm10 BigWigs (0.4-1 GB) at
10 kb bins: 69 s serial, 27 s parallel; 132k peak regions: 34 s → 12 s. The summary table is cached,
so changing only plot settings (method, colours, PCA vs correlation) reuses it in about a second.

Plot-specific:

| Key | Default | Plot | Notes |
|---|---|---|---|
| `method` | `pearson` | correlation | `spearman` |
| `corr_transform` | `raw` | correlation | `log1p` damps a few very high bins (Pearson) |
| `annotate` | true | correlation | print coefficients |
| `corr_cmap` | `RdBu_r` | correlation | |
| `pca_top` | 1000 | PCA | most variable features kept |
| `log` | `log2` | PCA | `raw` |
| `point_size` | 60 | PCA | |
| `maxrows` | 1000 | PairGrid | rows sampled for scatter, 10–50,000 |

PCA uses distinct markers per `group` only when groups pool samples, i.e. more than one group and
fewer groups than samples.

## Pitfalls

* Pearson on raw signal is dominated by a few huge bins. Prefer `spearman` or `corr_transform: log1p`
  for genome-wide bins.
* Small inputs (a single locus, or toy data) leave only a handful of bins, so the correlations mean
  little. Check the row count in `summary.tsv` and the "N retained" figure in `methods.txt`.
* A blacklist covering all the signal can leave fewer than 2 rows. The error then is "Fewer than two
  bins/regions have signal".
* `methods.txt` states the mode, bin size, chromosomes, skip list, blacklist and rows kept. Quote it.

Examples: `relationships-correlation-bins.json` (genome bins, no BED) and
`relationships-{correlation,pca,pairgrid}.json` (over peaks) in `skills/genomecanvas-plot/examples/`.

## Sample QC (large collections)

For tens to hundreds of samples with a sample sheet. It flags samples whose **mark** label or
**group** label does not fit their signal. The input is a precomputed bins × samples matrix, not
BigWigs. Unlike the other plots it reads no files from `files`; the sheet and matrix come from
settings.

Command line: `gcanvas-sampleqc` (installed by `pip install -e ".[studio]"`). Without the install,
use the module form, which is identical: `python3 -m genomecanvas.cli.sample_qc …`.

**Step 1: build the matrix (slow, run once, ideally on the cluster next to the BigWigs)**

```bash
gcanvas-sampleqc matrix --sheet samples.tsv --sample-col sample_id --path-col bw_path_hpc \
    --out nbl_10kb.npz --bin 10000 --workers 16 [--blacklist hg38-blacklist.bed] [--filter status=PASS]
```

* Each BigWig is read once. Per-sample checkpoints (`nbl_10kb.parts/`) make a killed job resume:
  rerun the same command. Missing or unreadable files are listed in `nbl_10kb.json` (exit code 1)
  and do not stop the run.
* Options: `--blacklist bed`, `--chroms chr1,chr2,…` (default: all but chrM, `*_random`, chrUn,
  alt/fix/hap/decoy), `--stat mean|mean0|sum|max`, `--checkpoint DIR`,
  `--path-map OLD=NEW` (rewrite a path prefix, e.g. `smb://server/groups=/Volumes/groups`).
* Cluster job (LSF; adapt queue and syntax to the scheduler):
  `bsub -n 16 -R "rusage[mem=2000]" -W 1:00 -o matrix.%J.log "gcanvas-sampleqc matrix … --workers 16"`.
* A deepTools `multiBigwigSummary bins|BED-file --outRawCounts` table works as the matrix too.

Measured cost, 4 hg38 BigWigs (300-515 MB) → 308,813 bins of 10 kb, Apple silicon, local SSD:

| Tool | Processes | Wall | CPU | Peak RSS (largest process) | Values |
|---|---|---|---|---|---|
| `gcanvas-sampleqc matrix` (pyBigWig, exact) | 4 | 17 s | 65 s | 445 MB | exact |
| deepTools `multiBigwigSummary` | 4 | 28-31 s | 101 s | 200-360 MB | zoom-level approximation (r ≥ 0.992 to exact) |
| bwtool `summary` | 4 | 53 s | 199 s | 37 MB | exact, text rounded to 2 decimals |

For 423 samples on 16 cores that is about 6-7 min for `gcanvas-sampleqc matrix` (expect 1.5-2x on a
network file system), about 11 min for deepTools and about 22 min for bwtool. The `.npz` is
0.2-0.4 GB at 10 kb.

**Step 2: run QC (seconds; repeat with other columns, filters and thresholds)**

```bash
gcanvas-sampleqc preview --sheet samples.tsv --matrix nbl_10kb.npz          # guesses, matches, tiers
gcanvas-sampleqc run --sheet samples.tsv --matrix nbl_10kb.npz --out qc/ \
    --mark target --group model --batch batch --filter "status=PASS"
```

`run` writes `qc_overview.png/pdf` instead of `figure.*`, plus the files listed under Outputs, and
prints status counts with the first 20 flagged samples as JSON. Only `run` needs repeating when
columns, filters or thresholds change.

From the runner, use `examples/relationships-sample-qc.json` (`qc_sheet`/`qc_matrix` are paths
locally, or file id/name/label on a server). In Studio, upload both files (they appear under
**Tables**), pick **Sample relationships → Sample QC**, choose them, then
**Preview & map columns…**. The preview shows the sheet, lets you pick the columns and filters,
reports how many rows match the matrix, and lists samples, groups, single-sample groups and the
size tier for each mark.

| Setting (Studio / runner) | CLI flag | Default | Notes |
|---|---|---|---|
| `qc_sheet`, `qc_matrix` | `--sheet`, `--matrix` | | required |
| `qc_sample` | `--sample-col` | best-matching column | its values must equal the matrix sample names (paths and `.bw` ignored) |
| `qc_mark`, `qc_group`, `qc_batch` | `--mark`, `--group`, `--batch` | guessed (`target`/`mark`, `group`/`model`/`condition`, `batch`) | batch optional |
| `qc_filters` | `--filter` (repeatable) | | `status=PASS; replicate!=NA`; `!=` excludes; several values comma-separated |
| `qc_skip_marks` | `--skip-marks` | `INPUT, IgG` | still mark-checked, not group-checked |
| `qc_k` | `--k` | 5 | neighbours voting on the mark |
| `qc_top` | `--top` | 5000 | most variable bins (per analysis) |
| `qc_method` | `--method` | `pearson` | on log(x+1); `spearman` |
| `qc_z` | `--z` | 3 | low-quality floor, robust SDs below the mark median |
| `qc_margin` | `--margin` | 0.02 | r by which another group must be closer before a sample is flagged |

How samples are judged:

1. **Mark check over all samples.** The k nearest samples by correlation vote. A sample is
   `mark_mismatch` when the vote names another mark and fewer than half of its neighbours share
   its own mark.
2. **Within each mark:**
   * `low_quality`: the sample's best r with anything **and** its best r with other groups both fall
     below the mark's median − z·MAD (Fisher z). These samples are then left out of everyone
     else's comparisons.
   * `mislabel`: its nearest member of another group beats its nearest own-group member by more
     than `qc_margin`. It becomes `swap` when a sample of that other group points back, and
     `mislabel_batch` when its nearest neighbour shares its batch (possibly a batch effect).
   * `undetermined`: no other usable member of its group in this mark, so it cannot be checked. It
     is never counted as passed.
   Nearest members are compared rather than group medians or centroids. With 1-3 samples per group,
   a centroid always looks closer than any single sample, and one swapped sample would drag its
   group-mates' medians down with it.

Outputs: `figure.png` (overview), `qc_flags.tsv` (one row per sample, most severe first, with
r values, nearest sample, closest other group and reasons), `qc_marks.tsv` (per mark: tier, groups,
single-sample groups, within/between-group median r, batch Δr), `qc_mark_confusion.tsv`,
`qc_group_matrix.tsv`, `qc_pca_by_mark.png` (one PCA per mark, flagged samples outlined) and
`qc_marks.pdf` (one page per mark: heatmap ordered by group with batch/group strips, group × group
median r, labelled PCA).

Reading the result:

* A swap between two groups with two members each implicates all four samples. The data cannot tell
  which pair moved, so check the IDs of all four.
* Marks where every group has one sample (NMYC in the NBL sheet; BRD4 and RNAPolII mostly)
  cannot be group-checked. Only the mark check applies there, and the warning names them.
* Flags are leads for a person to check against lab records, not verdicts. Report them as such.
