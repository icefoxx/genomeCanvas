# Heatmap pipeline (no deepTools)

Replacement for `deeptools_heatmaps.wrapper.py`. Matrices come from bwtool/pyBigWig;
figures are Matplotlib.

```bash
gcanvas heatmap \
  -p examples/heatmap_pipeline/oct4 \
  -m examples/heatmap_pipeline/meta.tsv \
  -r tests/fixtures/trackplot_demo/narrowpeak/H1_Oct4.bed,tests/fixtures/trackplot_demo/narrowpeak/H1_Nanog.bed \
  -l OCT4,NANOG \
  -f 1500 \
  --binsize 50 \
  --overwrite true
```

Python:

```python
from genomecanvas.pipeline import run_heatmap_pipeline
from genomecanvas.pipeline.run import HeatmapConfig
```
