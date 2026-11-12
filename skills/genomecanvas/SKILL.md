---
name: genomecanvas
description: Main entry for any GenomeCanvas request - turning ChIP-seq / ATAC / CUT&RUN signal (BigWig, bedGraph), peaks (BED, narrowPeak) and genes (GTF) into figures. Surveys the user's data, explains what is there and what is wrong with it, recommends figures with reasons, agrees a plan with the user, renders, looks at the result and refines it. Routes to the runner and the per-family skills. Start here whenever the user has genomics files and wants plots, QC or a figure for a paper.
---

# GenomeCanvas: from data to figures, in conversation

You are the user's analyst. Work in this loop, and talk to the user at the marked points:

```
goal → survey → [tell: what the data is] → recommend → [ask: agree the plan] → run → look → [show + refine] → report
```

All commands print one JSON object. Paths are relative to the repository root.

## 1. Goal

If the user named a question ("is SOX2 at enhancers?", "are my replicates OK?", "figure 2 for the
paper"), keep it in mind; it decides which recommendations matter. If not, do not interrogate them:
survey first and let the data suggest options.

## 2. Survey the data (always, before recommending anything)

```bash
python3 skills/genomecanvas-plot/scripts/gcplot.py survey DATA_DIR [MORE_FILES ...] --write PLAN_DIR
```

Reads headers, sample sheets and BED rows only, so it takes seconds even for many GB. It returns:

| Field | Use it to |
|---|---|
| `genome` | state the assembly (from BigWig chromosome sizes, else `project.yaml`) |
| `summary` | counts, marks/factors, controls (input/IgG), groups, total GB |
| `files[]` | per file: role, label, mark, replicate, group, control, GEO accession, region count, strongest peak |
| `skipped[]` | files it cannot use and why (e.g. `.gz` peaks must be decompressed) |
| `warnings[]` | problems: mixed assemblies, chr naming, >20,000 regions, no replicates, no control |
| `recommendations[]` | ranked figures: `plot`, `why`, `ask` (questions for the user), `cost`, `request` / `request_file` |

Labels and groups come from `samples.tsv` (columns `file label group`, optional `condition`,
`replicate`, `color`) and the genome / default locus from `project.yaml` when present; otherwise from
file names. If labels look like raw file names, offer to write a `samples.tsv`.

**Tell the user** in a few lines: assembly, what each sample is, controls and replicates, and every
warning in plain words. Name problems that limit the analysis (e.g. "one replicate per mark, so these
comparisons are descriptive"). Do not bury them.

## 3. Recommend and agree a plan

Offer 2-4 figures from `recommendations` (priority 1 first), each with one sentence of why, tied to
the user's goal when there is one. Put the recommendation's `ask` questions to the user, but only the
ones whose answer changes the figure; take sensible defaults for the rest and say which you took.

Ask, never guess:
* which locus/gene to show (suggest the survey's locus: project default or strongest peak);
* which groups form a contrast, and which is the reference;
* cluster names: render first, look, then propose names from what each block shows;
* whether to drop a control or a bad sample from a figure.

Decide yourself: bin sizes, flanks, colours, sort order, figure size (use defaults; change on request).

| The user wants … | Plot | Family skill |
|---|---|---|
| where marks/factors sit around peaks or TSS | `enrichment/profile-heatmap` (± `cluster`) | [enrichment](../genomecanvas-enrichment/SKILL.md) |
| region classes (promoter-like, enhancer-like) | `enrichment/profile-heatmap` with `cluster` | [enrichment](../genomecanvas-enrichment/SKILL.md) |
| sample similarity, replicate agreement, outliers | `relationships/correlation`, `pca` | [relationships](../genomecanvas-relationships/SKILL.md) |
| mislabelled / swapped samples in a large cohort | `relationships/sample-qc` | [relationships](../genomecanvas-relationships/SKILL.md#sample-qc-large-collections) |
| a browser-style view of one locus | `tracks/stacked` (+ GTF) | [tracks](../genomecanvas-tracks/SKILL.md) |
| ChIP over input, treatment vs control at a locus | `transform/log2-ratio`, `background` | [transform](../genomecanvas-transform/SKILL.md) |
| condition A vs B over peaks (≥2 replicates each) | `differential/volcano` (exploratory) | [differential](../genomecanvas-differential/SKILL.md) |
| one signal at two peak sets | `differential/peak-center` | [differential](../genomecanvas-differential/SKILL.md) |

Read the family skill before changing settings beyond the drafted request.

## 4. Run

```bash
python3 skills/genomecanvas-plot/scripts/gcplot.py run PLAN_DIR/01-....json --out results/01
python3 skills/genomecanvas-plot/scripts/gcplot.py run REQUEST.json --set clusters=5 --set group_names="A, B, C, D, E"
```

Run independent figures in parallel. Genome-wide tables take ~7 s per GB of BigWig the first time and
are cached afterwards; the survey's `cost` gives the estimate, so warn the user before long runs.
With a running Studio, add `--server URL --project ID` (see the
[runner skill](../genomecanvas-plot/SKILL.md) for request format, server mode and shared settings).

## 5. Look before you show

Open `figure.png` and check it as a reviewer would:

* Is anything empty, clipped or overlapping (legend over data, labels on thin clusters)?
* Does one sample (often input) or one huge cluster dominate? → drop the control, use
  `cluster_samples`, change `clusters` or `normalization`.
* Do the axes and colour ranges let the weaker samples be seen? (`zscale`, `y_range`)
* Do clusters mean something? Describe each block before naming it.

Fix obvious problems yourself and re-run; then show the figure. Say what it shows in one or two
sentences, in terms of the user's question, and what it cannot show.

## 6. Refine and report

Apply requested changes with `--set` on the same request (matrices are cached, so restyling takes
seconds). When done, give: the figure paths, `methods.txt` text for the paper, every `warnings`
entry, and the `request.json` that reproduces each figure.

## Function index

| Command | Purpose |
|---|---|
| `gcplot.py survey DIR --write PLAN` | classify data, flag problems, draft ranked requests |
| `gcplot.py list` | every plot slug with minimum inputs |
| `gcplot.py inspect FILE…` | role, chromosomes and extent of individual files |
| `gcplot.py run REQ [--set k=v] [--out DIR]` | render locally (in-process engine) |
| `gcplot.py run REQ --server URL --project ID` / `files --server …` | render through a running Studio |
| `gcanvas-sampleqc matrix / preview / run` | Sample QC for large collections ([details](../genomecanvas-relationships/SKILL.md#sample-qc-large-collections)) |
| `gcanvas-studio serve --config studio.toml` | run the multi-user web Studio |

Scripts live in `skills/genomecanvas-plot/scripts/`; `gcplot.py` = `python3 skills/genomecanvas-plot/scripts/gcplot.py`.
