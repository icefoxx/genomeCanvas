# Trackplot-parity figures (H1 ESC / OCT4)

These figures use the public demo files shipped with
[trackplot](https://github.com/PoisonAlien/trackplot): H1 embryonic-stem-cell
BigWigs and peaks around `chr6:31125776-31144789` (hg19, POU5F1/OCT4).

```bash
python -m examples.trackplot_parity.generate
```

| File | Contents |
|---|---|
| `locus.png` | Stacked OCT4, NANOG, H3K4me3, H3K4me1, H3K27ac, H2A.Z, CTCF, peaks, chromHMM, genes, ideogram |
| `locus_overlay.png` | Same signals as one overlaid line track |
| `profile.png` | Mean signal ±1.5 kb from OCT4/NANOG peak centers |
| `heatmap.png` | Region heatmap of the same windows |
| `pca.png` | Sample PCA on peak intensities |
| `volcano.png` | TF vs histone differential peak scores |
