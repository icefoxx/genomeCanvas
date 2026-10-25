# Architecture

The pipeline is raw genomic files → normalized data → tracks → FigureSpec → validation → layout →
backend-neutral rendering primitives → export. The current implementation includes deterministic
BedGraph/GTF/BED/VCF/GFF adapters, typed FigureSpec validation, signal transforms, Matplotlib vector
export, and provenance manifests. Optional BigWig querying is isolated behind pybigtools or pyBigWig, with bwtool used when present.
