# GenomeCanvas

GenomeCanvas is a Python engine for deterministic, reproducible genomic figures. It provides typed
figure specifications, genomic file readers, signal transformations, aligned track layouts,
Matplotlib rendering, and provenance manifests.

## Requirements

- Python 3.11 or newer
- `pyBigWig` for BigWig-based examples
- A Unix-like shell for the commands below; equivalent commands work on Windows PowerShell

## Build the Environment

Create an isolated virtual environment from the repository root:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,bigwig]'
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

The `bigwig` extra is optional for the core renderer, but is required by the real-data reference
plots. If BigWig support is not needed, install only `.[dev]`.

Check the installation:

```bash
python -c "import genomecanvas; print(genomecanvas.__version__)"
gcanvas --help
```

## Quickstart

Create a starter figure specification:

```bash
gcanvas init figure.yaml
```

Validate it:

```bash
gcanvas validate figure.yaml
```

Render a specification:

```bash
gcanvas render figure.yaml --output figure.pdf
```

The renderer writes the requested output and companion artifacts next to it:

- PDF, SVG, or PNG figure output
- `figure.spec.yaml`, the serialized figure specification
- `figure.manifest.json`, the provenance manifest

The generated starter specification has no tracks. To render a signal track, add a BedGraph source:

```yaml
version: '1.0'
genome:
  assembly: synthetic
region:
  chrom: chr1
  start: 0
  end: 300
tracks:
  - id: signal_a
    type: signal
    source:
      path: tests/fixtures/mini_genome/signal_a.bedgraph
```

## Figure Studio (web)

A browser studio for every plot type, from a single user up to a multi-user service for a
sequencing core (projects, accounts or SSO, admin console):

```bash
pip install -e ".[studio]"

# Personal: your own folder of BigWig/BED/GTF files, no sign-in → http://127.0.0.1:8765
python3 web/server.py --port 8765 --data-dir ~/data/bigwigs

# Multi-user: password accounts, projects, admin console at /admin
gcanvas-studio example-config > studio.toml        # set auth = "local", data_root, workspace
gcanvas-studio --config studio.toml user add you@core.edu --admin   # prints an invitation link
gcanvas-studio --config studio.toml project add P2026-017 --dir smith/P2026-017 --genome mm10
gcanvas-studio --config studio.toml project member P2026-017 client@uni.edu --role editor
gcanvas-studio --config studio.toml serve

# Production: Docker + nginx (+ optional SSO)
cd deploy && cp studio.toml studio.local.toml && docker compose up -d --build
```

Highlights: genome tracks, TSS/peak heatmaps and profiles, correlation/PCA over regions,
**genome-wide bins or sliding windows**, an interactive PCA (hover, draggable labels),
journal-size exports and an auto-generated Methods paragraph.

* [docs/STUDIO_QUICKSTART.md](docs/STUDIO_QUICKSTART.md): copy-paste recipes for the server and CLI
* [docs/STUDIO_DEPLOY.md](docs/STUDIO_DEPLOY.md): Docker/nginx/SSO, project manifests, sizing, security
* [web/README.md](web/README.md): features and HTTP API

## Generate p01-p30

The real-data fixtures are included in `tests/fixtures/reference_demo/`. Generate the main p01
figure and the p02-p30 comparison variants from the repository root:

```bash
python -m examples.reference_comparison.generate
python -m examples.reference_comparison.generate_p02_p30
```

The p01 PNG/PDF is written under `examples/reference_comparison/`. The p02-p30 variants are written
under `examples/reference_comparison/output/`. The scripts use the bundled GRCm38/mm10 BigWig,
GTF, and cytoband files, so no download is required.

Other runnable examples include:

```bash
python -m examples.demo_outputs.generate
python -m examples.advanced_features.generate
python -m examples.advanced_features.generate_comparison_plots
python -m examples.advanced_features.generate_real_data
```

## Development Checks

Run the test suite, lint, and compile checks with:

```bash
python -m pytest
ruff check src examples tests
python -m compileall -q src examples tests
```

Type checking is configured for the package:

```bash
mypy src
```

## Project Layout

```text
src/genomecanvas/       Core package, models, transforms, layouts, and renderers
examples/               Synthetic, advanced, and real-data figure generators
tests/                  Unit, integration, and bundled reference-data checks
docs/                   Architecture and scientific design documentation
```

Trackplot-parity extras (overlay, chromHMM, peaks, gene-symbol regions, profile
plots, heatmaps, PCA, and differential peaks) are documented in `docs/TRACK_API.md`.
Generate the mini demonstration with:

```bash
python -m examples.trackplot_parity.generate   # H1 ESC / OCT4 figures from trackplot demo data
```

The main architecture is:

```text
genomic files -> readers -> normalized tracks -> FigureSpec -> validation
              -> layout -> renderer -> exported figure and provenance
```

See `docs/ARCHITECTURE.md`, `docs/FIGURE_SPEC.md`, and `docs/TRACK_API.md` for implementation
details.
