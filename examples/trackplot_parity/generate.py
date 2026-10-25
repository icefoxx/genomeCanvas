"""Render trackplot-parity figures from the public H1 ESC demo data."""

from __future__ import annotations

from pathlib import Path

from genomecanvas.analysis import diffpeak, extract_summary, pca_plot, pca_scores, volcano_plot
from genomecanvas.model.figure_spec import FigureSpec
from genomecanvas.profile import extract_profile, profile_heatmap, profile_plot, profile_summarize
from genomecanvas.render.renderer import render

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/fixtures/trackplot_demo"
BW = DATA / "bw"
PEAKS = DATA / "narrowpeak"
OUT = Path(__file__).parent

# Same OCT4/POU5F1 window used in the trackplot README (hg19).
LOCUS = ("chr6", 31_125_776, 31_144_789)
SIGNALS = [
    ("oct4", "H1_Oct4.bw", "OCT4", "#d35400"),
    ("nanog", "H1_Nanog.bw", "NANOG", "#d35400"),
    ("k4me3", "H1_k4me3.bw", "H3K4me3", "#2980b9"),
    ("k4me1", "H1_k4me1.bw", "H3K4me1", "#2980b9"),
    ("k27ac", "H1_k27ac.bw", "H3K27ac", "#2980b9"),
    ("h2az", "H1_H2az.bw", "H2A.Z", "#27ae60"),
    ("ctcf", "H1_Ctcf.bw", "CTCF", "#27ae60"),
]


def locus_spec(*, overlay: bool = False) -> FigureSpec:
    """Build the H1 ESC OCT4 locus figure."""
    chrom, start, end = LOCUS
    tracks: list[dict[str, object]] = []
    member_ids: list[str] = []
    for track_id, filename, label, color in SIGNALS:
        member_ids.append(track_id)
        tracks.append({
            "id": track_id,
            "type": "signal",
            "label": label,
            "height": 0.9,
            "visibility": not overlay,
            "source": {"path": str(BW / filename)},
            "style": {"color": color, "ymin": 0},
        })
    if overlay:
        tracks.append({
            "id": "overlay",
            "type": "overlay",
            "label": "H1 ESC",
            "members": member_ids,
            "height": 2.2,
        })
    tracks.extend([
        {
            "id": "oct4_peaks",
            "type": "peak",
            "label": "OCT4 pk",
            "height": 0.5,
            "source": {"path": str(PEAKS / "H1_Oct4.bed")},
            "style": {"color": "#d35400"},
        },
        {
            "id": "nanog_peaks",
            "type": "peak",
            "label": "NANOG pk",
            "height": 0.5,
            "source": {"path": str(PEAKS / "H1_Nanog.bed")},
            "style": {"color": "#c0392b"},
        },
        {
            "id": "hmm",
            "type": "chromhmm",
            "label": "chromHMM",
            "height": 0.32,
            "source": {"path": str(PEAKS / "H1_chromHMM.bed")},
        },
        {
            "id": "genes",
            "type": "gene",
            "label": "genes",
            "height": 1.3,
            "style": {"genes": ["POU5F1", "TCF19"]},
        },
        {
            "id": "cyto",
            "type": "cytoband",
            "label": "chr6",
            "height": 0.4,
            "source": {"path": str(DATA / "cytoBandIdeo.txt.gz")},
        },
    ])
    return FigureSpec.parse_obj({
        "genome": {"assembly": "hg19"},
        "region": {"chrom": chrom, "start": start, "end": end},
        "gene_model": str(DATA / "chr6_refGene.gtf"),
        "collapse_transcripts": True,
        "highlights": [
            {"start": 31_139_576, "end": 31_140_109, "color": "#ffc41a", "alpha": 0.28},
            {"start": 31_140_574, "end": 31_141_112, "color": "#ffc41a", "alpha": 0.28},
        ],
        "scales": {"mode": "per_track"},
        "layout": {
            "width": 8.8,
            "height": 11.2 if not overlay else 7.4,
            "track_gap": 0.06,
            "show_y_axis": True,
            "y_tick_count": 3,
            "x_label": "hg19 chr6:31,125,776-31,144,789",
        },
        "tracks": tracks,
    })


def combined_peak_bed() -> Path:
    """Concatenate OCT4 and NANOG peaks for profile/heatmap/PCA windows."""
    path = OUT / "oct4_nanog_peaks.bed"
    path.write_text(
        (PEAKS / "H1_Oct4.bed").read_text(encoding="utf-8")
        + (PEAKS / "H1_Nanog.bed").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    return path


def main() -> None:
    if not (BW / "H1_Oct4.bw").is_file():
        raise SystemExit(f"missing trackplot demo BigWigs in {BW}")
    render(locus_spec(overlay=False), OUT / "locus.png")
    render(locus_spec(overlay=True), OUT / "locus_overlay.png")
    sources = {label: str(BW / filename) for _, filename, label, _ in SIGNALS}
    bed = combined_peak_bed()
    matrices = extract_profile(
        sources, bed=bed, start_from="center", up=1500, down=1500, binsize=50, input_type="bw",
    )
    colors = {label: color for _, _, label, color in SIGNALS}
    profile_plot(
        profile_summarize(matrices),
        up=1500,
        down=1500,
        output=OUT / "profile.png",
        colors=colors,
        xlabel="distance from peak center (bp)",
        ylabel="H1 ESC signal",
    )
    profile_heatmap(matrices, output=OUT / "heatmap.png", top_profile=True, zmax=8.0)
    summary = extract_summary(sources, bed, input_type="bw")
    pca_plot(pca_scores(summary, log2=True), output=OUT / "pca.png")
    conditions = {
        "OCT4": "TF",
        "NANOG": "TF",
        "CTCF": "TF",
        "H3K4me3": "histone",
        "H3K4me1": "histone",
        "H3K27ac": "histone",
        "H2A.Z": "histone",
    }
    volcano_plot(
        diffpeak(summary, conditions, numerator="TF", denominator="histone"),
        output=OUT / "volcano.png",
    )


if __name__ == "__main__":
    main()
