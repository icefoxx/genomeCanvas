# Quickstart

See the root `README.md` for the complete environment setup and usage instructions.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,bigwig]'
python -m pytest
python -m examples.reference_comparison.generate
python -m examples.reference_comparison.generate_p02_p30
```
