"""Sample QC: planted label errors are found, tables parse, the matrix builder is exact."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from genomecanvas.analysis.sample_qc import (
    SignalMatrix,
    apply_filters,
    guess_columns,
    load_matrix,
    match_samples,
    parse_filters,
    preview,
    qc_from_files,
    read_sheet,
    run_qc,
    write_tables,
)

FIX = Path(__file__).resolve().parent / "fixtures" / "trackplot_demo"


def simulate(seed=1, marks=4, groups=8, bins=12000):
    """Marks differ strongly, groups moderately (at signal bins), batches slightly."""
    rng = np.random.default_rng(seed)
    base = {f"M{m}": rng.gamma(0.6, 2.0, bins) for m in range(marks)}
    group_effect = {g: rng.normal(0, 0.6, bins) for g in range(groups)}
    batch_effect = {b: rng.normal(0, 0.15, bins) for b in range(3)}
    columns, meta = [], []
    for mark, signal in base.items():
        for g in range(groups):
            for r in range(1 + (g % 3)):  # 1, 2 or 3 replicates
                b = (g + r) % 3
                noise = rng.gamma(0.5, 0.4, bins)
                columns.append(
                    signal * np.exp(group_effect[g] * (signal > 1) + batch_effect[b]) + noise
                )
                meta.append(
                    dict(sample=f"{mark}_G{g}_r{r}", mark=mark, group=f"G{g}", batch=f"B{b}")
                )
    return np.column_stack(columns).astype(np.float32), meta


@pytest.fixture(scope="module")
def planted():
    values, meta = simulate()
    names = [m["sample"] for m in meta]
    index = {n: i for i, n in enumerate(names)}
    # Mark error: an M0 sample labelled M1.
    meta[index["M0_G1_r0"]]["mark"] = "M1"
    # Group swap inside M2 between two groups with 3 replicates (G2, G5).
    a, b = index["M2_G2_r0"], index["M2_G5_r0"]
    meta[a]["group"], meta[b]["group"] = "G5", "G2"
    # Noisy sample in M3, group G2 (3 replicates).
    c = index["M3_G2_r1"]
    values[:, c] = values[:, c] * 0.15 + np.random.default_rng(5).gamma(1.0, 2.0, values.shape[0])
    result = run_qc(SignalMatrix(names, values), meta)
    return result, {row["sample"]: row for row in result.rows}


def test_planted_errors_are_flagged(planted):
    result, rows = planted
    assert rows["M0_G1_r0"]["status"] == "mark_mismatch"
    assert rows["M0_G1_r0"]["predicted_mark"] == "M0"
    assert rows["M2_G2_r0"]["status"] == "swap" and rows["M2_G2_r0"]["swap_with"] == "M2_G5_r0"
    assert rows["M2_G5_r0"]["status"] == "swap"
    assert rows["M3_G2_r1"]["status"] == "low_quality"
    # The noisy replicate is left out of its group-mates' comparisons.
    assert rows["M3_G2_r0"]["status"] == "ok" and rows["M3_G2_r2"]["status"] == "ok"
    # Single-member groups are reported as not checkable, never as passing or failing.
    assert rows["M1_G0_r0"]["status"] == "undetermined"
    flagged = {r["sample"] for r in result.rows if r["status"] not in ("ok", "undetermined")}
    assert flagged == {"M0_G1_r0", "M2_G2_r0", "M2_G5_r0", "M3_G2_r1"}


def test_summaries_and_tables(planted, tmp_path):
    result, _ = planted
    marks = {m["mark"]: m for m in result.mark_summary}
    # G0, G3, G6 have one replicate; moving M0_G1_r0 to M1 leaves G1 with one too.
    assert marks["M0"]["singletons"] == 4 and marks["M0"]["tier"] == "S"
    assert result.confusion["M1"]["M0"] == 1
    names, grid, pairs = result.group_matrices["M1"]
    g0 = names.index("G0")
    assert np.isnan(grid[g0, g0]) and pairs[g0, g0] == 0  # one member: no within-group pair
    write_tables(result, tmp_path)
    flags = (tmp_path / "qc_flags.tsv").read_text().splitlines()
    assert flags[0].startswith("sample\tmark\tgroup\tbatch\tstatus")
    assert flags[1].split("\t")[4] == "mark_mismatch"  # most severe first
    assert {"qc_marks.tsv", "qc_mark_confusion.tsv", "qc_group_matrix.tsv"} <= {
        p.name for p in tmp_path.iterdir()
    }


def test_sheet_filters_guesses_and_matching(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\ttarget\tmodel\tbatch\tstatus\n"
        "A1\tH3K4me3\tX\tB1\tPASS\n"
        "A2\tH3K4me3\tY\tB1\tPASS\n"
        "A3\tH3K27ac\tX\tB2\tFAIL\n"
        "A4\tINPUT\tY\tB2\tPASS\n"
    )
    columns, rows = read_sheet(sheet)
    assert columns[:3] == ["sample_id", "target", "model"]
    assert len(apply_filters(rows, parse_filters("status=PASS; target!=INPUT"))) == 2
    with pytest.raises(ValueError, match="column=value"):
        parse_filters("status~PASS")
    with pytest.raises(ValueError, match="not in the sample sheet"):
        apply_filters(rows, parse_filters("nope=1"))
    samples = ["/data/A2.bw", "'A1.bigWig'", "A4"]
    assert match_samples(rows, "sample_id", samples) == {0: 1, 1: 0, 3: 2}
    guesses = guess_columns(columns, rows, samples)
    assert guesses == dict(sample="sample_id", mark="target", group="model", batch="batch")


def test_load_matrix_formats(tmp_path):
    raw = tmp_path / "counts.tab"
    raw.write_text(
        "#'chr'\t'start'\t'end'\t'A1.bw'\t'A2.bw'\nchr1\t0\t10\t1\tnan\nchr1\t10\t20\t3\t4\n"
    )
    m = load_matrix(raw)
    assert m.samples == ["A1.bw", "A2.bw"] and m.values.tolist() == [[1, 0], [3, 4]]
    np.savez(tmp_path / "m.npz", values=np.ones((3, 2)), samples=np.array(["a", "b"]), bin_size=10)
    assert load_matrix(tmp_path / "m.npz").bin_size == 10


def test_qc_from_files_and_preview(tmp_path):
    values, meta = simulate(seed=2, marks=2, groups=4, bins=4000)
    sheet = tmp_path / "s.tsv"
    sheet.write_text(
        "sample_id\ttarget\tmodel\tbatch\tstatus\n"
        + "".join(f"{m['sample']}\t{m['mark']}\t{m['group']}\t{m['batch']}\tPASS\n" for m in meta)
        + "extra\tM0\tG0\tB0\tFAIL\n"
    )
    matrix = tmp_path / "m.npz"
    np.savez(matrix, values=values, samples=np.array([m["sample"] + ".bw" for m in meta]))
    out = preview(sheet, matrix, {"qc_filters": "status=PASS"})
    assert out["chosen"]["mark"] == "target" and out["n_matched"] == len(meta)
    assert {m["mark"] for m in out["marks"]} == {"M0", "M1"}
    result, notes = qc_from_files(sheet, matrix, {"qc_filters": "status=PASS"})
    assert notes["matched"] == len(meta) and not notes["unmatched"]
    assert Counter(r["status"] for r in result.rows)["mark_mismatch"] == 0
    with pytest.raises(ValueError, match="match matrix columns"):
        qc_from_files(sheet, matrix, {"qc_sample": "model"})


def test_matrix_builder_matches_pybigwig_and_resumes(tmp_path):
    pyBigWig = pytest.importorskip("pyBigWig")
    from genomecanvas.pipeline.sample_matrix import build_matrix

    sources = {n: str(FIX / "bw" / f"H1_{n}.bw") for n in ("k4me3", "Ctcf")}
    sources["ghost"] = str(tmp_path / "missing.bw")
    out = tmp_path / "m.npz"
    report = build_matrix(sources, out, size=1000, chroms=["chr6"], workers=2, log=lambda _: None)
    assert report["missing_files"] == {"ghost": sources["ghost"]} and report["samples"] == 2
    m = load_matrix(out)
    rows = np.flatnonzero((m.starts >= 31125000) & (m.starts < 31145000))
    handle = pyBigWig.open(sources["Ctcf"])
    expected = handle.stats("chr6", 31125000, 31145000, type="mean", nBins=20, exact=True)
    assert np.allclose(
        m.values[rows, m.samples.index("Ctcf")], np.nan_to_num(np.array(expected, float))
    )
    logs = []
    build_matrix(sources, out, size=1000, chroms=["chr6"], workers=2, log=logs.append)
    assert all("checkpoint" in line for line in logs)
