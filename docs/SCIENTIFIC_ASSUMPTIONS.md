# Scientific assumptions

- Internal coordinates are always 0-based, half-open intervals.
- GTF/GFF/VCF coordinates are converted explicitly from 1-based closed conventions.
- Signal differences require matching regions and declared units; ratio operations require an explicit pseudocount.
- No biological normalization or differential testing is inferred by the renderer.
- BedGraph bins are display aggregations and never modify source files.
