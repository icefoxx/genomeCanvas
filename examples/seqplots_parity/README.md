# SeqPlots example parity

Uses the public SeqPlots `inst/extdata` files (ce10 chrI 100 kb).

```bash
python examples/seqplots_parity/generate.py
```

Outputs land in `output/`: locus, transcript heatmap/profile, peak heatmap/profile, and a two-feature-set heatmap matching `getPlotSetArray(bw, c(transcripts, peaks), "ce10")`.
