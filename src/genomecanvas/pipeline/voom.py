"""Extract Up/Down/Control BEDs from a differential-peak VOOM TSV."""

from __future__ import annotations

from pathlib import Path

THRESHOLDS = {
    "Up2FC_Down2FC": ("p.value", 1.1, 1.0, "Up2FC", "Down2FC"),
    "Up2_Down2": ("q.value", 0.05, 1.0, "Up2", "Down2"),
    "Up2NoFDR_Down2NoFDR": ("p.value", 0.05, 1.0, "Up2NoFDR", "Down2NoFDR"),
    "Up_Down": ("q.value", 0.05, 0.0, "Up", "Down"),
    "UpNoFDR_DownNoFDR": ("p.value", 0.05, 0.0, "UpNoFDR", "DownNoFDR"),
    "Control": ("p.value", 1.0, 0.0, "NA", "NA"),
}


def extract_voom_beds(
    path: str | Path,
    prefix: str | Path,
    *,
    threshold: str = "Up2_Down2",
    include_control: bool = False,
    seed: int = 7201,
) -> list[tuple[str, str]]:
    """Write BED files for the requested regulation classes."""
    key = threshold.replace("_Control", "")
    include_control = include_control or threshold.endswith("_Control") or threshold == "Control"
    if key not in THRESHOLDS and threshold != "Control":
        raise ValueError(f"unsupported voom threshold {threshold!r}")
    rows = _read_voom(path)
    stat_field, stat_cut, fc_cut, up_name, down_name = THRESHOLDS.get(key, THRESHOLDS["Up_Down"])
    up_rows = [
        row for row in rows
        if row.get(stat_field, 1.0) < stat_cut and row.get("log2FC", 0.0) > fc_cut
    ] if "Up" in threshold else []
    down_rows = [
        row for row in rows
        if row.get(stat_field, 1.0) < stat_cut and row.get("log2FC", 0.0) < -fc_cut
    ] if "Down" in threshold else []
    ctrl_rows = [row for row in rows if row.get("Regulation") == "Control"]
    beds: list[tuple[str, str]] = []
    base = Path(prefix)
    base.parent.mkdir(parents=True, exist_ok=True)
    if len(up_rows) > 10:
        dest = Path(f"{base}.{up_name}.bed")
        _write_bed(dest, up_rows)
        beds.append((str(dest), f"{up_name} [{len(up_rows)} regions]"))
    if len(down_rows) > 10:
        dest = Path(f"{base}.{down_name}.bed")
        _write_bed(dest, down_rows)
        beds.append((str(dest), f"{down_name} [{len(down_rows)} regions]"))
    if include_control and ctrl_rows:
        n_ctrl = min(len(ctrl_rows), max(len(up_rows), len(down_rows), 1))
        rng = np_random(seed)
        chosen = [ctrl_rows[i] for i in rng.choice(len(ctrl_rows), n_ctrl, replace=False)]
        dest = Path(f"{base}.Control.bed")
        _write_bed(dest, chosen)
        beds.append((str(dest), f"Control [{len(chosen)} regions]"))
    return beds


def _read_voom(path: str | Path) -> list[dict[str, float | str]]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    header = [item.strip() for item in lines[0].split("\t")]
    rows: list[dict[str, float | str]] = []
    for raw in lines[1:]:
        if not raw.strip():
            continue
        fields = raw.split("\t")
        record = {
            key: fields[index] if index < len(fields) else ""
            for index, key in enumerate(header)
        }
        region = str(record.get("Region", ""))
        chrom, rest = region.split(":") if ":" in region else (record.get("chrm", ""), "0-1")
        start_s, end_s = rest.split("-") if "-" in rest else ("0", "1")
        parsed: dict[str, float | str] = {
            "chrom": chrom,
            "start": int(start_s.replace(",", "")),
            "end": int(end_s.replace(",", "")),
            "Regulation": str(record.get("Regulation", "")),
            "log2FC": _to_float(record.get("log2FC", 0)),
            "p.value": _to_float(record.get("p.value", 1)),
            "q.value": _to_float(record.get("q.value", 1)),
        }
        rows.append(parsed)
    return rows


def _write_bed(path: Path, rows: list[dict[str, float | str]]) -> None:
    lines = [f"{row['chrom']}\t{int(row['start'])}\t{int(row['end'])}" for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _to_float(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0.0


def np_random(seed: int):
    import numpy as np

    return np.random.default_rng(seed)
