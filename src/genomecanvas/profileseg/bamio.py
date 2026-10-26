"""BAM input for ProfileSeg (v3): fragment counts per bin around anchors, with read-level filters.

Requires pysam (``pip install pysam``; optional dependency ``genomecanvas[bam]``).

Signal modes
  midpoint   one count per fragment at its midpoint (MNase-seq dyads, CUT&RUN, ATAC)
  coverage   fragment coverage (ChIP-seq style), as base-pair overlap per bin
  end5/end3  one count per read at its 5′ or 3′ end (PRO-seq, run-on, exo assays), strand-aware

Paired-end fragments are counted once (from the leftmost proper-pair mate, using TLEN). Single-end
reads are extended to ``fragment_length``; when it is not given it is estimated from the strand
cross-correlation of read 5′ ends over the anchors (the distance between the + and − strand peaks;
cf. FunChIP's compute_fragments_length).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

MODES = ("midpoint", "coverage", "end5", "end3")


@dataclass
class BamOptions:
    mode: str = "midpoint"
    min_mapq: int = 10
    proper_pairs: bool = True
    remove_duplicates: bool = True
    min_fragment: int = 0
    max_fragment: int = 1000
    fragment_length: int | None = None  # single-end extension; None → estimated
    strand: str = "any"  # end modes: any | sense | antisense (relative to the anchor strand)


@dataclass
class BamSample:
    counts: np.ndarray  # regions × bins (float; integer for midpoint/end modes)
    fragments_per_region: np.ndarray  # fragments with a midpoint inside the window
    library_size: float  # fragments (or reads) passing filters on the indexed contigs
    paired: bool
    fragment_length: float | None
    fragment_hist: np.ndarray = field(default_factory=lambda: np.zeros(0))  # 0..max_fragment, paired-end only


def _require_pysam():
    try:
        import pysam  # noqa: F401
    except ImportError as err:  # pragma: no cover - depends on the environment
        raise ImportError("BAM input needs pysam: pip install pysam (or genomecanvas[bam])") from err
    import pysam
    return pysam


def _keep(read, opt: BamOptions) -> bool:
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    if read.mapping_quality < opt.min_mapq:
        return False
    if opt.remove_duplicates and read.is_duplicate:
        return False
    return True


def detect_paired(bam, n: int = 2000) -> bool:
    k = paired = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        paired += read.is_paired
        k += 1
        if k >= n:
            break
    return k > 0 and paired / k > 0.5


def estimate_fragment_length(bam, windows: list[tuple[str, int, int]], opt: BamOptions,
                             lo: int = 50, hi: int = 400) -> float:
    """Single-end fragment length from the strand cross-correlation of read 5′ ends."""
    plus, minus = [], []
    for chrom, start, end in windows[:2000]:
        if chrom not in bam.references:
            continue
        n = end - start
        p, m = np.zeros(n), np.zeros(n)
        for read in bam.fetch(chrom, max(start, 0), end):
            if not _keep(read, opt):
                continue
            if read.is_reverse:
                x = read.reference_end - 1 - start
                if 0 <= x < n:
                    m[x] += 1
            else:
                x = read.reference_start - start
                if 0 <= x < n:
                    p[x] += 1
        plus.append(p); minus.append(m)
    if not plus:
        return float(opt.fragment_length or 200)
    best, best_d = -np.inf, 200
    for d in range(lo, hi + 1, 5):
        score = sum(float(np.dot(p[: len(p) - d], m[d:])) for p, m in zip(plus, minus, strict=True))
        if score > best:
            best, best_d = score, d
    return float(best_d)


def library_size(bam, opt: BamOptions, paired: bool) -> float:
    """Mapped reads on indexed contigs (halved for paired-end): a library-size proxy for CPM scaling."""
    total = sum(s.mapped for s in bam.get_index_statistics())
    return total / 2.0 if paired else float(total)


def read_bam(path: str, regions: list[dict], anchors: np.ndarray, edges: np.ndarray, oriented: bool,
             opt: BamOptions | None = None) -> BamSample:
    """Counts per bin over anchor + edges for every region (minus-strand windows mirrored when oriented)."""
    pysam = _require_pysam()
    opt = opt or BamOptions()
    if opt.mode not in MODES:
        raise ValueError(f"bam mode must be one of {MODES}")
    bam = pysam.AlignmentFile(path, "rb")
    paired = detect_paired(bam)
    lo_e, hi_e = float(edges[0]), float(edges[-1])
    pad = opt.max_fragment if paired else int(opt.fragment_length or 400)
    windows = [(r["chrom"], int(a + lo_e) - pad, int(a + hi_e) + pad) for r, a in zip(regions, anchors, strict=True)]
    frag_len = None
    if not paired and opt.mode in ("midpoint", "coverage"):
        frag_len = float(opt.fragment_length) if opt.fragment_length else estimate_fragment_length(bam, windows, opt)
    nb = len(edges) - 1
    counts = np.zeros((len(regions), nb))
    frags = np.zeros(len(regions))
    hist = np.zeros(opt.max_fragment + 1)
    seen_hist = set()
    for k, (r, a) in enumerate(zip(regions, anchors, strict=True)):
        chrom = r["chrom"]
        if chrom not in bam.references:
            continue
        minus = oriented and r.get("strand") == "-"
        rel_list, cov_list = [], []
        for read in bam.fetch(chrom, max(windows[k][1], 0), max(windows[k][2], 1)):
            if not _keep(read, opt):
                continue
            if opt.mode in ("end5", "end3"):
                if opt.strand != "any":
                    same = read.is_reverse == (r.get("strand") == "-")
                    if (opt.strand == "sense") != same:
                        continue
                five = read.reference_end - 1 if read.is_reverse else read.reference_start
                three = read.reference_start if read.is_reverse else read.reference_end - 1
                rel_list.append((five if opt.mode == "end5" else three) - a)
                continue
            if paired:
                if opt.proper_pairs and not read.is_proper_pair:
                    continue
                tlen = read.template_length
                if tlen <= 0:  # count each fragment once, from its leftmost mate
                    continue
                if tlen < opt.min_fragment or tlen > opt.max_fragment:
                    continue
                start, end = read.reference_start, read.reference_start + tlen
                if read.query_name not in seen_hist:  # a fragment near two anchors is counted once
                    hist[tlen] += 1
                    seen_hist.add(read.query_name)
            else:
                if read.is_reverse:
                    end = read.reference_end
                    start = end - int(frag_len)
                else:
                    start = read.reference_start
                    end = start + int(frag_len)
            if opt.mode == "midpoint":
                rel_list.append((start + end) // 2 - a)
            else:
                cov_list.append((start - a, end - a))
        if rel_list:
            rel = np.asarray(rel_list, float)
            if minus:
                rel = -rel
            inside = (rel >= lo_e) & (rel < hi_e)
            counts[k] = np.histogram(rel[inside], bins=edges)[0]
            frags[k] = inside.sum()
        if cov_list:
            iv = np.asarray(cov_list, float)
            if minus:
                iv = -iv[:, ::-1]
            # overlap of each fragment with each bin (bp), divided by the bin width
            ov = np.clip(np.minimum(iv[:, 1:2], edges[1:]) - np.maximum(iv[:, 0:1], edges[:-1]), 0, None)
            counts[k] = ov.sum(0) / np.diff(edges)
            mids = iv.mean(1)
            frags[k] = np.sum((mids >= lo_e) & (mids < hi_e))
    lib = library_size(bam, opt, paired)
    bam.close()
    return BamSample(counts=counts, fragments_per_region=frags, library_size=lib, paired=paired,
                     fragment_length=frag_len, fragment_hist=hist if paired else np.zeros(0))


def fragment_length_summary(hist: np.ndarray) -> dict:
    """Median, interquartile range and share of nucleosome-sized (120–180 bp) fragments."""
    if hist.sum() == 0:
        return {}
    x = np.arange(len(hist))
    cdf = np.cumsum(hist) / hist.sum()
    q = lambda p: float(x[np.searchsorted(cdf, p)])  # noqa: E731
    return dict(median=q(0.5), q25=q(0.25), q75=q(0.75),
                nucleosomal_share=float(hist[120:181].sum() / hist.sum()))
