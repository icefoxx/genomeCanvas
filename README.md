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

The main architecture is:

```text
genomic files -> readers -> normalized tracks -> FigureSpec -> validation
              -> layout -> renderer -> exported figure and provenance
```

See `docs/ARCHITECTURE.md`, `docs/FIGURE_SPEC.md`, and `docs/TRACK_API.md` for implementation
details.
