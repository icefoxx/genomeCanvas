# Local background-adjusted visualization

GenomeCanvas supports a visualization-only local background transformation.
The transformation estimates background from explicit flanking regions and
materializes a `DerivedTrack` containing the parent track, estimate, operation,
pseudocount, and result signal.

```yaml
transformations:
  - id: local_enrichment
    operation: log2_ratio
    input: WT_ATAC
    background:
      strategy: flanking
      upstream: 10000
      downstream: 10000
      estimator: median
    pseudocount: 0.01
```

Supported operations are `subtract`, `ratio`, and `log2_ratio`. The median is
the default estimator because it is less sensitive to a nearby peak. Mean and
trimmed mean are explicit alternatives.

This is **not sample normalization**. It does not alter BigWig files, infer
biological normalization, or support differential testing. Use labels such as
“local background-adjusted signal” or “log2(local enrichment)”, not “normalized
signal”.
