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

Feature: Trackplot-parity locus, profile, and analysis plots
Public behavior observed: IGV-style locus tracks, chromHMM bars, peak overlays, gene-symbol regions, profile plots, region heatmaps, PCA, and two-group peak tests are common epigenomic figure types.
Independent design decision: Extend typed FigureSpec track types and Python APIs; query BigWig through pybigtools/pyBigWig with optional bwtool; keep display transforms and tests independent of the R implementation.
Implementation files: `src/genomecanvas/io`, `src/genomecanvas/profile`, `src/genomecanvas/analysis`, `src/genomecanvas/render`
Tests: `tests/test_trackplot_parity.py`
