# Clean-room development log

This project is independently designed and implemented. No source code, tests, documentation, or
distinctive internal design was copied from genomic plotting engines.

Future entries use this format:

```text
Feature:
Public behavior observed:
Independent design decision:
Implementation files:
Tests:
```

Feature: Declarative FigureSpec and deterministic genomic tracks
Public behavior observed: Structured locus figures commonly combine signal, annotations, genes, and comparisons.
Independent design decision: Use typed Pydantic models, normalized interval objects, backend-neutral primitives, and explicit transformations.
Implementation files: `src/genomecanvas/model`, `src/genomecanvas/tracks`, `src/genomecanvas/render`
Tests: `tests/test_scaffold.py`, `tests/test_features.py`
