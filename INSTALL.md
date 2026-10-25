# Installation

GenomeCanvas requires Python 3.11 or newer. For the complete setup, see the root `README.md`.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,bigwig]'
```

Use `.[dev]` instead when BigWig support is not required.

The `bigwig` extra installs pyBigWig and pybigtools. Put `bwtool` on PATH to use it as a
BigWig engine. Select the backend with `engine="bwtool"` / `engine="pybigwig"`, YAML
`bigwig_engine`, `gcanvas --bigwig-engine`, or `GENOMECANVAS_BIGWIG_ENGINE`.
