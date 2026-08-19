# Advanced extensions

GenomeCanvas now includes the TASK-101–108 extension layer:

- deterministic derived transformations and parent metadata;
- flanking and rolling local-background estimation;
- subtract, ratio, log2-ratio, and local z-score operations;
- MAD, robust quantile, shared robust, and symmetric display scaling;
- scientific figure lint diagnostics;
- automatic largest-peak and enrichment focus selection;
- priority label filtering with optional leader-line instructions;
- manifest transformation history.

All derived operations are visualization transformations. They do not modify
source signal files and must not be described as sample normalization or
differential testing.
