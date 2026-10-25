# GenomeCanvas Studio: web front end

Browser studio for all GenomeCanvas plot types, backed by the Python engine in
`src/genomecanvas/studio/`. This folder holds the static UI: `index.html`, `app.js` and `styles.css`
for the studio, `login.*` for sign-in and invitations, and `admin.*` for the admin console.

## Run

Single user, local data (no sign-in):

```sh
python3 web/server.py --port 8765 --data-dir ~/data/my_bigwigs
```

Multi-user for clients of a core facility (projects, accounts or SSO, admin console):

```sh
gcanvas-studio --config studio.toml serve
```

See **[docs/STUDIO_DEPLOY.md](../docs/STUDIO_DEPLOY.md)** for Docker/nginx/SSO deployment,
project manifests (`samples.tsv`, `project.yaml`), administration and sizing.

## What clients can do

* **Projects**: choose a project in the top bar. Files delivered by the core arrive pre-labelled
  (labels, groups and colours from `samples.tsv`). Extra files can be dropped onto the page and
  count against the project's quota.
* **Plot types**:
  * Genomic tracks: stacked or overlaid tracks, gene models, chromHMM states, highlights.
  * Enrichment: heatmap, profile, profile + heatmap, clustered heatmap; anchored at the
    centre, TSS or TES, or scale-regions. Minus-strand regions are oriented like deepTools.
  * Sample relationships: correlation, PCA, PairGrid over **selected regions**,
    **genome-wide bins** or **sliding windows**, with blacklist, chromosome filters,
    statistic, zero filtering and subsampling.
  * Differential analysis (exploratory) and signal transformations.
* **Interactive PCA**: hover a point for sample, group, PC1–3 and file; drag labels to tidy them.
  Static exports use automatic label repulsion. A standalone interactive HTML can be downloaded.
* **Styling**: 10 curated palettes or a custom one; colour-scale galleries; fonts; axis styles;
  exact export size in mm or inches (journal column presets); PNG up to 1200 dpi; editable text in
  SVG/PDF; an *Advanced* JSON panel for any documented style field.
* **Reproducibility**: every result includes `methods.txt` (a draft Methods paragraph),
  `provenance.json` (input hashes, software versions) and `request.json`. **Save** keeps a figure
  together with its settings in the project for collaborators. Templates apply a core-standard
  figure in one click.
* **Productivity**: ⌘/Ctrl+Enter to plot, **Live** mode, `/` to search settings, gene search and
  locus navigation, a history strip of recent figures, and configuration import/export.

## HTTP API (summary)

All state-changing requests need the header `X-GenomeCanvas: studio` and same-origin.

| method & path | purpose |
|---|---|
| `POST /api/login`, `/api/logout`, `/api/invite` | password accounts |
| `GET /api/me`, `/api/meta` | current user, projects, catalogue/palettes |
| `GET /api/projects/<id>` | project, files, templates |
| `POST /api/projects/<id>/uploads` | raw body, URL-encoded `X-Filename` |
| `DELETE /api/projects/<id>/files/<file>` | delete an upload |
| `GET /api/projects/<id>/genes?files=…&q=…` | gene search |
| `POST /api/projects/<id>/jobs`, `GET /api/jobs/<job>` | render a figure, poll status |
| `GET /results/<job>/<file>[?download=name]` | outputs (access-checked) |
| `GET/POST /api/projects/<id>/figures`, `PATCH/DELETE /api/figures/<id>` | saved figures |
| `POST /api/templates`, `DELETE /api/templates/<id>` | templates (admin) |
| `/api/admin/…` | users, projects, members, scans, audit, cleanup (admin) |
| `GET /vendor/plotly.min.js` | Plotly bundle from the installed `plotly` package |
