# Installation

GenomeCanvas requires Python 3.11 or newer. For the complete setup, see the root `README.md`.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,bigwig]'
```

Use `.[dev]` instead when BigWig support is not required.
