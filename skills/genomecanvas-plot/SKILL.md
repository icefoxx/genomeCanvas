---
name: genomecanvas-plot
description: Render any GenomeCanvas Studio figure (genomic tracks, enrichment heatmaps/profiles, correlation/PCA/PairGrid, differential volcano/peak-center, signal transformations) from BigWig/bedGraph, BED and GTF files, either locally in-process or by submitting a job to a running Studio web server (agent mode). Runner and request reference; for an open-ended request ("plot my data") start with the genomecanvas skill, which surveys the data and plans the figures.
---

# GenomeCanvas plotting (runner)

One runner, `scripts/gcplot.py`, drives every Studio plot. It takes the same request the web UI
sends, applies the same UI defaults, and writes the same outputs (`figure.png/svg/pdf`,
`methods.txt`, `provenance.json`, `request.json`, plus plot-specific tables). Every command prints a
single JSON object on stdout; warnings from libraries go to stderr.

Paths below are relative to the repository root.

## 1. Pick the plot

| Plot slug | Minimum inputs | Skill with the parameters |
|---|---|---|
| `tracks/stacked`, `tracks/overlay` | 1 / 2 signals | [genomecanvas-tracks](../genomecanvas-tracks/SKILL.md) |
| `tracks/genes` | 1 GTF | [genomecanvas-tracks](../genomecanvas-tracks/SKILL.md) |
| `enrichment/heatmap`, `profile`, `profile-heatmap`, `clustered-heatmap` | 1 signal + 1 BED | [genomecanvas-enrichment](../genomecanvas-enrichment/SKILL.md) |
| `relationships/correlation`, `pca`, `pairgrid` | 2 signals (+ 1 BED only in `regions` mode) | [genomecanvas-relationships](../genomecanvas-relationships/SKILL.md) |
| `relationships/sample-qc` | sample sheet + bins × samples matrix (no `files`) | [genomecanvas-relationships](../genomecanvas-relationships/SKILL.md#sample-qc-large-collections) |
| `differential/volcano` | 4 signals (2 groups × ≥2 replicates) + 1 BED | [genomecanvas-differential](../genomecanvas-differential/SKILL.md) |
| `differential/peak-center` | 1 signal + 2 BED | [genomecanvas-differential](../genomecanvas-differential/SKILL.md) |
| `transform/difference`, `ratio`, `log2-ratio`, `background` | 2 signals | [genomecanvas-transform](../genomecanvas-transform/SKILL.md) |

`python3 skills/genomecanvas-plot/scripts/gcplot.py list` prints this table from the engine itself.
`gcplot.py survey DATA_DIR --write PLAN_DIR` classifies a folder and drafts ranked requests; the
[genomecanvas](../genomecanvas/SKILL.md) skill describes the full survey → plan → run → review loop.
Read the family skill before choosing settings: several defaults and input rules are non-obvious.
Sample QC also has its own command line (`gcanvas-sampleqc matrix | preview | run`) for building
the matrix on a cluster; see the relationships skill.

## 2. Write a request

```json
{
  "plot": "relationships/correlation",
  "title": "H1 marks, 10 kb bins",
  "files": [
    {"path": "data/H1_k4me3.bw", "role": "signal", "label": "H3K4me3", "group": "Histone"},
    {"path": "data/H1_k27ac.bw", "role": "signal", "label": "H3K27ac", "group": "Histone"}
  ],
  "settings": {"rel_mode": "bins", "rel_bin": 10000, "method": "spearman"},
  "advanced": {}
}
```

* `files[].role` is `signal` (`.bw .bigwig .bedgraph .bdg .bg`), `regions` (`.bed .narrowPeak
  .broadPeak`), `genes` (`.gtf .gff`) or `table` (`.tsv .csv .npz`); it is guessed from the suffix
  when omitted.
* `label` is what the figure shows (default: file name stem). `group` drives PCA markers, the
  volcano contrast and heatmap group strips. `color` (hex) and `cmap` override the palette per file.
* Order matters: palette colours follow file order, and transforms use the first two signals.
* Relative paths resolve against the request file's folder.
* Ready-made requests for every plot on the bundled demo data are in `examples/`.

## 3. Run it

Local (in-process engine, local paths):

```bash
python3 skills/genomecanvas-plot/scripts/gcplot.py run request.json --out results/fig1
python3 skills/genomecanvas-plot/scripts/gcplot.py run request.json --set dpi=600 --set theme=minimal
python3 skills/genomecanvas-plot/scripts/gcplot.py inspect data/*.bw peaks.bed   # roles, chromosomes, extent
```

Agent / server mode (same request, submitted to a running Studio):

```bash
python3 skills/genomecanvas-plot/scripts/gcplot.py files --server http://127.0.0.1:8765 --project demo
python3 skills/genomecanvas-plot/scripts/gcplot.py run request.json --server http://127.0.0.1:8765 --project demo
```

* In server mode each `files[]` entry is matched to a project file by `id`, then file `name` (or the
  basename of `path`), then `label`. `rel_blacklist` is resolved the same way.
* Project defaults (`project.yaml` / demo defaults such as the locus) are applied underneath the
  request settings, as when the UI opens the project.
* Sign-in deployments need the session cookie: `--cookie <gc_session value>` or env
  `GCANVAS_SESSION`. The runner sends the required `X-GenomeCanvas: studio` and `Origin` headers.
* The server enforces its own limits (48 files per request, per-user job limit, queue size,
  project access). An agent running inside the web server should go through this jobs API rather
  than import the engine, so authorisation, quotas and worker isolation still apply.

Result JSON (both modes): `status`, `out`, `figure`, `files`, `methods`, `warnings`, `counts`,
`elapsed`; server mode adds `job` and `result_urls`. On failure: `{"status": "error", "error": ...}`
and exit code 1. Engine errors are written for end users; relay them as they are.

## 4. Settings shared by every plot

| Key | Default | Notes |
|---|---|---|
| `title` | plot name | |
| `palette` | `studio` | `studio okabe npg lancet nejm tol morandi nord tableau ink` |
| `custom_palette` | | `"#264653, #2a9d8f"`; overrides `palette` |
| `theme` | `clean` | `clean classic grid minimal` |
| `font_family` | DejaVu Sans | any installed font; SVG/PDF text stays editable |
| `title_size`, `axis_label_size`, `tick_size` | 12, 10, 8 | points |
| `size_unit` | `mm` | `mm` or `in` |
| `fig_width`, `fig_height` | auto | journal columns: 85 / 114 / 180 mm |
| `exact_size` | false | export exactly width × height (needs both) |
| `dpi` | 300 | PNG only, 72–1200 |
| `transparent` | false | |
| `linewidth`, `fill_alpha` | per family | |

`advanced` accepts only `rc` (Matplotlib rcParams under `axes. xtick. ytick. grid. lines. font.
legend. …`), `heatmap`, `profile`, `layout`, `style`, `track_style`. Anything else is rejected.

## 5. Report back

Show the figure (`figure.png`), quote `methods.txt` when the user will write up the figure, and pass
on every entry in `warnings`. They carry caveats the user needs, e.g. that the volcano test is exploratory.
Keep `request.json` with the figure; re-running it reproduces the figure.
