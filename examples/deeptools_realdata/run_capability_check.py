"""Run GenomeCanvas on real deepTools / GEO BigWig + BED inputs."""

from __future__ import annotations

import traceback
from pathlib import Path

from genomecanvas.analysis import extract_summary, pca_plot, pca_scores
from genomecanvas.io.bigwig import bigwig_chroms
from genomecanvas.io.signal import query_signal
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.pipeline.bins import summarize_bins, write_genomic_bins
from genomecanvas.pipeline.cluster import cluster_labels
from genomecanvas.pipeline.correlation import correlation_heatmap
from genomecanvas.pipeline.matrix import compute_region_matrix, sort_matrix
from genomecanvas.pipeline.metadata import SampleMeta
from genomecanvas.pipeline.plot_style import HeatmapPlotStyle, ProfilePlotStyle
from genomecanvas.pipeline.plotting import plot_enrichment_heatmap, plot_profile_lines
from genomecanvas.render.renderer import render
from genomecanvas.tracks.difference import compare_signal

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).parent / "output"
PGT = ROOT / "tests/fixtures/deeptools_encode"
OFFICIAL = PGT / "official"
H1 = ROOT / "tests/fixtures/trackplot_demo"
GEO = ROOT / "tests/fixtures/reference_demo"
OUT.mkdir(exist_ok=True)

REPORT: list[str] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    mark = "YES" if ok else "NO "
    line = f"{mark}  {name}" + (f" — {detail}" if detail else "")
    REPORT.append(line)
    print(line, flush=True)


def main() -> None:
    fly_bws = [
        PGT / "bw/pgt_chrx_2e6_5e6.bw",
        PGT / "bw/pgt_chrx_2.5e6_3.5e6.bw",
        PGT / "bw/pgt_chrx_2.5e6_3.5e6_b.bw",
    ]
    windows = PGT / "bed/chrX_windows_200.bed"
    fly_genes = PGT / "bed/dm3_chrX_2_5Mb_refGene.bed"

    # 1) Official deepTools computeMatrix inputs
    try:
        samples = [
            SampleMeta(str(OFFICIAL / "computeMatrix2.bw"), "dt_example", "rep", "chip", "Blues", "auto")
        ]
        matrix = compute_region_matrix(
            samples, [(str(OFFICIAL / "computeMatrix2.bed"), "genes")],
            flank=50, binsize=10, ref_point="center", skip_zeros=False, engine="auto",
        )
        plot_enrichment_heatmap(
            matrix, OUT / "official_heatmap.png",
            style=HeatmapPlotStyle(
                title="deepTools computeMatrix2 example",
                labels=("coverage",),
                xlab="distance from center",
                ylab="features",
                colorbar_label="signal",
                height=5.0,
                width_per_sample=3.2,
            ),
        )
        plot_profile_lines(
            matrix, OUT / "official_profile.png",
            style=ProfilePlotStyle(title="Mean profile", labels=("coverage",), xlab="distance from center"),
        )
        record(
            "Official deepTools computeMatrix2.bw + BED heatmap",
            True,
            f"{matrix.values.shape[0]} regions × {matrix.values.shape[1]} bins",
        )
    except Exception as exc:
        record("Official deepTools computeMatrix2.bw + BED heatmap", False, str(exc))
        traceback.print_exc()

    # 2) Locus IGV-style on pyGenomeTracks chrX window
    try:
        tracks = []
        colors = ["#388b83", "#bb8063", "#5087bc"]
        for i, path in enumerate(fly_bws):
            tracks.append({
                "id": path.stem, "type": "signal", "label": path.stem,
                "source": {"path": str(path)}, "style": {"color": colors[i], "ymin": 0},
            })
        spec = FigureSpec.parse_obj({
            "genome": {"assembly": "dm3"},
            "region": {"chrom": "chrX", "start": 2_500_000, "end": 3_500_000},
            "tracks": tracks,
            "layout": {"width": 9, "height": 5.5, "x_label": "chrX", "x_tick_unit": "Mb"},
            "style": {"title": "pyGenomeTracks chrX BigWigs (deepTools org test data)"},
        })
        render(spec, OUT / "fly_locus.png")
        record("Locus tracks from real BigWig", True, "chrX:2.5-3.5 Mb, 3 bw")
    except Exception as exc:
        record("Locus tracks from real BigWig", False, f"{type(exc).__name__}: {exc}")
        traceback.print_exc()

    # 3) computeMatrix-style reference-point heatmap (200 windows)
    try:
        samples = [
            SampleMeta(str(path), path.stem, "rep", "chip", "Blues" if i else "Oranges", "auto")
            for i, path in enumerate(fly_bws)
        ]
        matrix = compute_region_matrix(
            samples, [(str(windows), "windows")],
            flank=500, binsize=50, ref_point="center", skip_zeros=True, engine="auto",
        )
        matrix = sort_matrix(matrix, method="mean")
        plot_enrichment_heatmap(
            matrix, OUT / "fly_heatmap.png",
            style=HeatmapPlotStyle(
                title="chrX windows · pyGenomeTracks",
                labels=("2–5 Mb", "2.5–3.5 Mb A", "2.5–3.5 Mb B"),
                xlab="distance from window center",
                ylab="windows",
                colorbar_label="coverage",
                height=7.0,
            ),
        )
        plot_profile_lines(
            matrix, OUT / "fly_profile.png",
            style=ProfilePlotStyle(
                title="Mean coverage around window centers",
                labels=("2–5 Mb", "2.5–3.5 Mb A", "2.5–3.5 Mb B"),
                xlab="distance from window center",
                ylab="mean coverage",
            ),
        )
        record(
            "Region heatmap + profile (computeMatrix reference-point)",
            True,
            f"{matrix.values.shape[0]} regions × {len(samples)} samples",
        )
    except Exception as exc:
        record("Region heatmap + profile (computeMatrix reference-point)", False, str(exc))
        traceback.print_exc()
        matrix = None

    # 4) kmeans clustered heatmap
    try:
        if matrix is None:
            raise RuntimeError("no matrix")
        labels = cluster_labels(matrix, method="kmeans", n_clusters=3, seed=7)
        record("K-means clustering on enrichment matrix", True, f"clusters={sorted(set(int(x) for x in labels))}")
    except Exception as exc:
        record("K-means clustering on enrichment matrix", False, str(exc))

    # 5) correlation / PCA (multiBigwigSummary BED-file analogue)
    try:
        summary = extract_summary(
            {p.stem: str(p) for p in fly_bws}, windows, engine="auto"
        )
        correlation_heatmap(summary.values, summary.samples, OUT / "fly_corr.png")
        pca_plot(pca_scores(summary, log2=True), output=OUT / "fly_pca.png")
        record(
            "Correlation + PCA from region means",
            True,
            f"{summary.values.shape}, finite={bool(summary.values.size and (summary.values == summary.values).all())}",
        )
    except Exception as exc:
        record("Correlation + PCA from region means", False, str(exc))
        traceback.print_exc()

    # 6) scale-regions on matching dm3 genes
    try:
        samples = [
            SampleMeta(str(fly_bws[0]), "pgt", "rep", "chip", "Blues", "auto")
        ]
        scaled = compute_region_matrix(
            samples, [(str(fly_genes), "genes")],
            flank=500, binsize=50, ref_point="TSS", body_length=1000,
            skip_zeros=True, engine="auto",
        )
        plot_enrichment_heatmap(
            scaled, OUT / "scale_regions.png",
            style=HeatmapPlotStyle(
                title="dm3 refGene scale-regions (TSS–TES)",
                labels=("chrX coverage",),
                xlab="TSS → TES",
                ylab="genes",
                colorbar_label="coverage",
                height=8.0,
            ),
        )
        record(
            "Scale-regions (TSS + body) on dm3 refGene BED",
            True,
            f"shape={scaled.values.shape}",
        )
    except Exception as exc:
        record("Scale-regions (TSS + body) on dm3 refGene BED", False, f"{type(exc).__name__}: {exc}")
        traceback.print_exc()

    # 7) difference of two real BigWigs (chrX vs X aliases)
    try:
        a = query_signal(str(fly_bws[0]), "chrX", 2_500_000, 3_500_000, pixels=800)
        b = query_signal(str(fly_bws[1]), "X", 2_500_000, 3_500_000, pixels=800)
        diff = compare_signal(a, b, operation="difference")
        record("BigWig difference (bamCompare analogue on tracks)", True, f"n={diff.values.size}")
    except Exception as exc:
        record("BigWig difference (bamCompare analogue on tracks)", False, str(exc))
        traceback.print_exc()

    # 8) Official bamCoverage phiX locus + difference
    try:
        left = OFFICIAL / "bamCoverage_result4.bw"
        right = OFFICIAL / "bamCoverage_result5.bw"
        spec = FigureSpec.parse_obj({
            "genome": {"assembly": "phiX174"},
            "region": {"chrom": "phiX174", "start": 0, "end": 5386},
            "tracks": [
                {
                    "id": "cov4", "type": "signal", "label": "bamCoverage 4",
                    "source": {"path": str(left)}, "style": {"color": "#1f77b4", "ymin": 0},
                },
                {
                    "id": "cov5", "type": "signal", "label": "bamCoverage 5",
                    "source": {"path": str(right)}, "style": {"color": "#d62728", "ymin": 0},
                },
            ],
            "layout": {"width": 8, "height": 3.5},
            "style": {"title": "Official deepTools bamCoverage BigWigs (phiX174)"},
        })
        render(spec, OUT / "official_bamcoverage_locus.png")
        a = query_signal(str(left), "phiX174", 0, 5386, pixels=200)
        b = query_signal(str(right), "phiX174", 0, 5386, pixels=200)
        compare_signal(a, b, operation="log2_ratio", pseudocount=0.1)
        record("Official bamCoverage BigWig locus + log2 ratio", True, "phiX174")
    except Exception as exc:
        record("Official bamCoverage BigWig locus + log2 ratio", False, str(exc))
        traceback.print_exc()

    # 9) Human H1 ESC (trackplot) matching peaks
    try:
        h1 = [
            H1 / "bw/H1_k4me3.bw",
            H1 / "bw/H1_k27ac.bw",
            H1 / "bw/H1_Oct4.bw",
        ]
        samples = [
            SampleMeta(str(p), p.stem, "H1", "h1", "Blues", "auto") for p in h1
        ]
        beds = [
            (str(H1 / "narrowpeak/H1_Oct4.bed"), "H1_OCT4"),
            (str(H1 / "narrowpeak/H1_k4me3.bed"), "H1_H3K4me3"),
        ]
        human = compute_region_matrix(
            samples, beds, flank=1000, binsize=50, ref_point="center",
            skip_zeros=False, engine="auto",
        )
        human = sort_matrix(human, method="mean")
        plot_enrichment_heatmap(
            human, OUT / "human_h1_peaks.png",
            style=HeatmapPlotStyle(
                title="H1 ESC ChIP at OCT4 / H3K4me3 peaks",
                labels=("H3K4me3", "H3K27ac", "OCT4"),
                xlab="distance from peak center",
                ylab="peaks",
                colorbar_label="signal",
                height=6.0,
            ),
        )
        record(
            "Heatmap: H1 ESC BigWig × matching H1 peak BEDs",
            True,
            f"{human.values.shape[0]} peaks",
        )
    except Exception as exc:
        record("Heatmap: H1 ESC BigWig × matching H1 peak BEDs", False, str(exc))
        traceback.print_exc()

    # 10) Human locus with genes
    try:
        spec = FigureSpec.parse_obj({
            "genome": {"assembly": "hg19"},
            "region": {"chrom": "chr6", "start": 31_125_776, "end": 31_144_789},
            "gene_model": str(H1 / "chr6_refGene.gtf"),
            "collapse_transcripts": True,
            "tracks": [
                {
                    "id": "k4me3", "type": "signal", "label": "H3K4me3",
                    "source": {"path": str(H1 / "bw/H1_k4me3.bw")},
                    "style": {"color": "#2980b9", "ymin": 0},
                },
                {
                    "id": "k27ac", "type": "signal", "label": "H3K27ac",
                    "source": {"path": str(H1 / "bw/H1_k27ac.bw")},
                    "style": {"color": "#27ae60", "ymin": 0},
                },
                {
                    "id": "oct4", "type": "signal", "label": "OCT4",
                    "source": {"path": str(H1 / "bw/H1_Oct4.bw")},
                    "style": {"color": "#8e44ad", "ymin": 0},
                },
                {
                    "id": "peaks", "type": "peak",
                    "source": {"path": str(H1 / "narrowpeak/H1_Oct4.bed")},
                    "style": {"color": "#d35400"},
                },
                {"id": "genes", "type": "gene"},
            ],
            "layout": {"width": 9, "height": 6, "x_tick_unit": "auto", "transcript_label": "name"},
            "style": {"title": "H1 ESC ChIP-seq at OCT4/POU5F1 (trackplot demo)"},
        })
        render(spec, OUT / "human_locus.png")
        record("Human locus (H1 BigWig + peaks + RefSeq genes)", True)
    except Exception as exc:
        record("Human locus (H1 BigWig + peaks + RefSeq genes)", False, str(exc))
        traceback.print_exc()

    # 11) GEO mouse 2-cell + HoxD regulatory BED
    try:
        geo_bws = [
            GEO / "GSM5746910_MS_2cell_Input_rep1.bigWig",
            GEO / "GSM5746912_MS_2cell_IP_rep1.bigWig",
        ]
        if all(path.is_file() for path in geo_bws):
            samples = [
                SampleMeta(str(geo_bws[0]), "Input_rep1", "input", "input", "Greys", "auto"),
                SampleMeta(str(geo_bws[1]), "IP_rep1", "ip", "chip", "Reds", "auto"),
            ]
            hoxd = PGT / "bed/HoxD_cluster_regulatory_regions_mm10.bed"
            mouse = compute_region_matrix(
                samples, [(str(hoxd), "HoxD")],
                flank=2000, binsize=100, ref_point="center",
                skip_zeros=False, engine="auto",
            )
            mouse = sort_matrix(mouse, method="mean")
            plot_enrichment_heatmap(
                mouse, OUT / "mouse_hoxd_heatmap.png",
                style=HeatmapPlotStyle(
                    title="GEO GSM5746910/12 at HoxD regulatory regions",
                    labels=("2-cell Input", "2-cell IP"),
                    xlab="distance from region center",
                    ylab="HoxD regions",
                    colorbar_label="signal",
                    height=6.0,
                ),
            )
            spec = FigureSpec.parse_obj({
                "genome": {"assembly": "mm10"},
                "region": {"chrom": "chr2", "start": 74_668_309, "end": 74_765_142},
                "gene_model": str(PGT / "bed/HoxD.gtf"),
                "collapse_transcripts": True,
                "tracks": [
                    {
                        "id": "input", "type": "signal", "label": "2-cell Input",
                        "source": {"path": str(geo_bws[0])},
                        "style": {"color": "#7f8c8d", "ymin": 0},
                    },
                    {
                        "id": "ip", "type": "signal", "label": "2-cell IP",
                        "source": {"path": str(geo_bws[1])},
                        "style": {"color": "#c0392b", "ymin": 0},
                    },
                    {
                        "id": "reg", "type": "peak",
                        "source": {"path": str(hoxd)},
                        "style": {"color": "#f39c12"},
                    },
                    {"id": "genes", "type": "gene"},
                ],
                "layout": {"width": 10, "height": 5.5, "x_tick_unit": "auto", "transcript_label": "name"},
                "style": {"title": "GEO GSM5746910/12 at HoxD cluster (mm10)"},
            })
            render(spec, OUT / "mouse_hoxd_locus.png")
            record(
                "GEO mouse 2-cell BigWig × HoxD BED heatmap + locus",
                True,
                f"{mouse.values.shape[0]} regions",
            )
        else:
            record("GEO mouse 2-cell BigWig × HoxD BED heatmap + locus", False, "GEO BigWigs missing")
    except Exception as exc:
        record("GEO mouse 2-cell BigWig × HoxD BED heatmap + locus", False, str(exc))
        traceback.print_exc()

    # 12) multiBigwigSummary bins analogue
    try:
        chroms = bigwig_chroms(str(fly_bws[0]))
        bin_bed = write_genomic_bins(chroms, 200_000, OUT / "fly_bins_200kb.bed")
        samples = [
            SampleMeta(str(path), path.stem, "rep", "chip", "Blues", "auto")
            for path in fly_bws
        ]
        labels, table = summarize_bins(samples, bin_bed, engine="auto")
        correlation_heatmap(table, labels, OUT / "fly_bin_corr.png")
        record(
            "multiBigwigSummary bins (genome-wide tiles)",
            True,
            f"{table.shape[0]} bins × {table.shape[1]} samples",
        )
    except Exception as exc:
        record("multiBigwigSummary bins (genome-wide tiles)", False, str(exc))
        traceback.print_exc()

    record("bamCoverage / bamCompare from BAM", False, "engine has no BAM reader (by design)")
    record("plotFingerprint / plotCoverage QC", False, "not implemented")

    (OUT / "CAPABILITY.md").write_text(
        "# GenomeCanvas vs real deepTools-style inputs\n\n"
        + "\n".join(f"- {r}" for r in REPORT)
        + "\n",
        encoding="utf-8",
    )
    print("\nWrote", OUT / "CAPABILITY.md")


if __name__ == "__main__":
    main()
