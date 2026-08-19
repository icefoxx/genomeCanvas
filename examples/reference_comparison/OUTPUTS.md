# p02-p30 real-data outputs

Generate all outputs with:

```bash
python -m examples.reference_comparison.generate_p02_p30
```

Each reference-style name has an independently rendered `.png` and `.pdf` in
this directory's `output/` folder. The output families use the downloaded
BigWig/GRCm38 inputs, actual cytobands, and real transcript annotations. They
are behaviorally comparable rather than pixel copies of the reference images.

`p10_ratio_08` and `p11_ratio_15` use y/x figure aspect ratios of 0.8 and 1.5.
`p12_yscale_gene` uses one y-range per locus across tracks; `p13_yscale_track`
uses one y-range per track across loci.
