"""ProfileSeg: metric properties, extraction, inference design, acceptance scenarios, per-region tests."""

from __future__ import annotations

import gzip
import json

import numpy as np
import pytest

from genomecanvas.profileseg import metrics as M
from genomecanvas.profileseg.analysis import analyse, exact_permutation, holm, permutation_assignments
from genomecanvas.profileseg.benchmark import run_benchmark
from genomecanvas.profileseg.data import Config, Profiles, load_config, read_deeptools_matrix
from genomecanvas.profileseg.regional import bh, half_sample_mode, moderated_t, per_region_tests
from genomecanvas.profileseg.simulate import SimSettings, simulate

EDGES = np.arange(-2000, 2001, 20.0)
C = M.centers(EDGES)


def gauss(mu, sigma, c=C):
    return np.exp(-0.5 * ((c - mu) / sigma) ** 2)


# ---------------------------------------------------------------- metric properties (plan unit tests)


def test_w1_identity_translation_rescaling_and_w2():
    p = M.normalise(gauss(-100, 120))
    q = M.normalise(gauss(100, 120))
    assert M.wasserstein1(p, p, EDGES) == pytest.approx(0, abs=1e-12)
    assert M.wasserstein1(p, q, EDGES) == pytest.approx(200, rel=1e-6)
    assert M.wasserstein2(p, q, EDGES) == pytest.approx(200, rel=1e-3)
    assert M.wasserstein1(
        M.normalise(5 * gauss(0, 150)), M.normalise(gauss(0, 150)), EDGES
    ) == pytest.approx(0, abs=1e-9)
    assert 0 <= M.jensen_shannon(p, q) <= 1 and M.jensen_shannon(p, p) == pytest.approx(
        0, abs=1e-12
    )


def test_w1_with_nonuniform_bins_uses_physical_coordinates():
    edges = np.concatenate([np.arange(-2000, 0, 40.0), np.arange(0, 2001, 10.0)])
    c, w = M.centers(edges), M.widths(edges)
    p = M.normalise(gauss(-300, 100, c) * w)
    q = M.normalise(gauss(-100, 100, c) * w)
    assert M.wasserstein1(p, q, edges) == pytest.approx(200, abs=0.5)


def test_asymmetry_reverses_under_orientation_reversal_and_moments():
    p = M.normalise(gauss(150, 100) + 0.3 * gauss(-300, 80))
    assert M.asymmetry(p[::-1], EDGES) == pytest.approx(-M.asymmetry(p, EDGES))
    g = M.normalise(gauss(0, 150))
    assert M.centroid(g, EDGES) == pytest.approx(0, abs=1e-9)
    assert M.spread_sd(g, EDGES) == pytest.approx(150, rel=0.01)
    assert M.iqr_width(g, EDGES) == pytest.approx(1.349 * 150, rel=0.01)


def test_background_estimators_and_amplitude_units():
    values = np.full((1, 1, len(C)), 2.0)
    values[0, 0, 90:110] = 10.0
    assert M.background(values, EDGES, 500)[0, 0] == pytest.approx(2.0)
    amp = M.amplitude_table(values, EDGES, center_half=250, flank_width=500)
    assert amp["auc"][0, 0] == pytest.approx(2.0 * 4000 + 8.0 * 20 * 20)  # signal × bp
    neg = values - 5
    assert M.amplitude_table(neg, EDGES, 250, 500)["negative_fraction"][0, 0] > 0.5


# ---------------------------------------------------------------- inference design


def test_permutation_resolution_matches_design():
    labels, _ = permutation_assignments(3, 3)
    assert len(labels) == 20 and labels[0].tolist() == [False] * 3 + [True] * 3
    _, floor = exact_permutation(np.array([1.0, 1.1, 0.9, 3.0, 3.1, 2.9]), labels)
    assert floor == pytest.approx(0.1)  # two-sided 3 vs 3: 0.05 is unattainable
    labels2, _ = permutation_assignments(2, 2)
    assert exact_permutation(np.array([1.0, 1.2, 5.0, 5.2]), labels2)[1] == pytest.approx(1 / 3)
    paired, _ = permutation_assignments(3, 3, pairs=[(0, 3), (1, 4), (2, 5)])
    assert len(paired) == 8


def test_inclusion_is_condition_independent():
    profiles, _ = simulate("amplitude", SimSettings(n_regions=120), seed=3)
    a = analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=20))
    swapped = Profiles(
        profiles.values,
        profiles.edges,
        [{**s, "condition": "KO" if s["condition"] == "WT" else "WT"} for s in profiles.samples],
        profiles.regions,
        profiles.oriented,
    )
    b = analyse(swapped, Config(reference="WT", alternative="KO", bootstrap=20))
    assert np.array_equal(a.included, b.included)


def test_unoriented_anchors_drop_signed_metrics():
    profiles, _ = simulate("shift", SimSettings(n_regions=80, oriented=False), seed=4)
    a = analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=20))
    names = {e["metric"] for e in a.effects}
    assert "centroid" not in names and "asymmetry" not in names and "sd_width" in names


def test_degenerate_inputs_fail_clearly():
    zero = Profiles(
        np.zeros((4, 10, len(C))),
        EDGES,
        [dict(sample_id=f"s{i}", condition="WT" if i < 2 else "KO") for i in range(4)],
        [dict(chrom="c", start=i, end=i + 1, name=f"r{i}", strand="+") for i in range(10)],
        True,
    )
    with pytest.raises(ValueError, match="regions pass inclusion"):
        analyse(zero, Config(reference="WT", alternative="KO", bootstrap=10))
    with pytest.raises(ValueError, match="need samples"):
        analyse(zero, Config(reference="WT", alternative="MUT", bootstrap=10))


def test_missing_bins_exclude_region_not_crash():
    profiles, _ = simulate("null", SimSettings(n_regions=60), seed=5)
    profiles.values[2, 7, :30] = np.nan  # e.g. window beyond a chromosome end in one sample
    a = analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=20))
    assert not a.included[7] and a.exclusion["missing_bins"] == 1


# ---------------------------------------------------------------- acceptance: each case described correctly


@pytest.fixture(scope="module")
def bench():
    return run_benchmark(
        scenarios=["null", "amplitude", "shift", "broaden", "right_flank", "left_flank"],
        reps=8,
        settings=SimSettings(n_regions=200),
        bootstrap=50,
        seed=4242,
    )


def _rate(bench, scenario, metric):
    return next(r for r in bench["summary"] if r["scenario"] == scenario and r["metric"] == metric)


def test_amplitude_only_is_not_reported_as_shape_change(bench):
    assert _rate(bench, "amplitude", "log2_auc")["welch_flag_rate"] >= 0.9
    for m in ("sd_width", "iqr_width", "centroid"):
        assert _rate(bench, "amplitude", m)["welch_flag_rate"] <= 0.25, m
    assert (
        _rate(bench, "amplitude", "naive_sd_width")["welch_flag_rate"] >= 0.9
    )  # the confound, exposed


def test_shift_width_and_flank_cases(bench):
    shift = _rate(bench, "shift", "centroid")
    assert shift["welch_flag_rate"] >= 0.9 and shift["direction_accuracy"] >= 0.9
    assert shift["mean_estimate"] == pytest.approx(100, abs=8)
    assert _rate(bench, "shift", "log2_auc")["welch_flag_rate"] <= 0.25
    assert _rate(bench, "broaden", "iqr_width")["welch_flag_rate"] >= 0.9
    assert _rate(bench, "broaden", "centroid")["welch_flag_rate"] <= 0.25
    assert _rate(bench, "right_flank", "asymmetry")["mean_estimate"] > 0.1
    assert _rate(bench, "left_flank", "asymmetry")["mean_estimate"] < -0.1


def test_null_false_positive_rate_is_bounded(bench):
    for m in ("log2_auc", "centroid", "sd_width", "iqr_width", "asymmetry"):
        assert _rate(bench, "null", m)["welch_flag_rate"] <= 0.25, m


# ---------------------------------------------------------------- inputs and end-to-end


def test_config_sections_and_unknown_keys(tmp_path):
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "contrast: {reference: A, alternative: B}\nbackground: {width: 400, estimator: median}\n"
        "shape: {core_flank: 800}\nflank: 2000\nbinsize: 20\n"
    )
    c = load_config(cfg)
    assert (
        c.reference,
        c.alternative,
        c.background_width,
        c.background_estimator,
        c.core_flank,
    ) == ("A", "B", 400, "median", 800)
    cfg.write_text("flank: 2000\nbogus: 1\n")
    with pytest.raises(ValueError, match="Unknown ProfileSeg setting"):
        load_config(cfg)


def test_deeptools_matrix_reader(tmp_path):
    header = {
        "upstream": [100],
        "downstream": [100],
        "body": [0],
        "bin size": [20],
        "ref point": ["TSS"],
        "sample_labels": ["WT_1", "KO_1"],
        "sample_boundaries": [0, 10, 20],
    }
    path = tmp_path / "m.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@" + json.dumps(header) + "\n")
        handle.write("chr1\t100\t200\tg1\t0\t+\t" + "\t".join(["1"] * 10 + ["2"] * 10) + "\n")
        handle.write("chr1\t500\t600\tg2\t0\t-\t" + "\t".join(["3"] * 10 + ["nan"] * 10) + "\n")
    p = read_deeptools_matrix(
        path, [dict(sample_id="WT_1", condition="WT"), dict(sample_id="KO_1", condition="KO")]
    )
    assert p.values.shape == (2, 2, 10) and p.oriented and np.isnan(p.values[1, 1]).all()
    assert p.edges[0] == -100 and p.edges[-1] == 100


def test_example_end_to_end_with_bigwig_extraction(tmp_path):
    pytest.importorskip("pyBigWig")
    from genomecanvas.profileseg.cli import main

    out = tmp_path / "ex"
    assert main(["example", "--out", str(out), "--n-regions", "60", "--bootstrap", "50"]) == 0
    results = out / "results"
    for rel in (
        "reports/00_discovery.md",
        "reports/01_methods_and_assumptions.md",
        "reports/02_QC.md",
        "reports/04_limitations.md",
        "results/metrics/effects.tsv",
        "reports/05_regional.md",
        "results/metrics/region_tests.tsv",
        "figures/E_effect_sizes.pdf",
        "figures/G_region_volcano.pdf",
        "run_log.json",
    ):
        assert (results / rel).is_file(), rel
    effects = (results / "results/metrics/effects.tsv").read_text()
    centroid = next(line for line in effects.splitlines() if line.startswith("centroid"))
    assert 70 < float(centroid.split("\t")[3]) < 130  # shift scenario: +100 bp
    assert (
        "smallest attainable two-sided p is 0.1"
        in (results / "reports/01_methods_and_assumptions.md").read_text()
    )


def test_extraction_orients_and_mirrors_minus_strand(tmp_path):
    pyBigWig = pytest.importorskip("pyBigWig")
    from genomecanvas.profileseg.data import extract

    path = tmp_path / "s.bw"
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([("chr1", 20000)])
    starts = list(range(8000, 12000, 20))
    bw.addEntries(
        ["chr1"] * len(starts),
        starts,
        ends=[s + 20 for s in starts],
        values=[float(i) for i in range(len(starts))],
    )  # ramp: increasing to the right
    bw.close()
    regions = [
        dict(chrom="chr1", start=10000, end=10001, name="p", strand="+"),
        dict(chrom="chr1", start=10000, end=10001, name="m", strand="-"),
    ]
    samples = [dict(sample_id="a", condition="WT", signal_file=str(path))]
    cfg = Config(
        anchor="start",
        orientation="auto",
        flank=1000,
        binsize=20,
        background_width=200,
        core_flank=600,
    )
    prof = extract(cfg, samples, regions)
    plus, minus = prof.values[0, 0], prof.values[0, 1]
    assert prof.oriented and np.all(np.diff(plus) > 0) and np.all(np.diff(minus) < 0)
    assert np.allclose(
        plus, minus[::-1] - (minus[::-1] - plus).mean()
    )  # mirror image of the same ramp


# ---------------------------------------------------------------- v1.1: roles, multiplicity, diagnostics


def test_holm_and_roles():
    adj = holm(np.array([0.01, 0.04, 0.03, np.nan]))
    assert adj[:3] == pytest.approx([0.03, 0.06, 0.06]) and np.isnan(adj[3])
    profiles, _ = simulate("global_scale", SimSettings(n_regions=200), seed=21)
    a = analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=20))
    roles = {e["metric"]: e["role"] for e in a.effects}
    assert roles["log2_auc"] == "primary" and roles["naive_sd_width"] == "comparator"
    calls = {e["metric"]: e["call"] for e in a.effects}
    assert calls["log2_auc"] == "higher" and calls["log2_enrichment"].startswith("no clear change")
    assert any("local background differs" in f for f in a.flags)  # global scaling is visible in QC


def test_no_spurious_flags_without_problems():
    for scenario in ("null", "amplitude", "shift"):
        a = analyse(simulate(scenario, SimSettings(n_regions=200), seed=31)[0],
                    Config(reference="WT", alternative="KO", bootstrap=20))
        assert not a.flags, (scenario, a.flags)


def test_opposite_shifts_are_flagged_as_possibly_bidirectional():
    a = analyse(simulate("opposite_shift", SimSettings(), seed=7)[0],
                Config(reference="WT", alternative="KO", bootstrap=20))
    assert any("bidirectional" in f for f in a.flags)


# ---------------------------------------------------------------- v2: per-region moderated tests


def test_moderated_t_matches_limma():
    # Reference values from limma 3.58.1: eBayes(lmFit(t(x), model.matrix(~g))) on the same matrix.
    i, j = np.arange(1, 7)[:, None], np.arange(1, 41)[None, :]
    x = np.sin(i * j * 0.7) + 0.3 * np.cos(j * 1.3) * (i > 3) + 0.05 * j * (i % 2)
    r = moderated_t(x, np.arange(3), np.arange(3, 6))
    assert r["d0"] == pytest.approx(8.4257293121, rel=1e-8)
    assert r["s02"] == pytest.approx(0.8968050434, rel=1e-8)
    assert r["t"][[0, 6, 22]] == pytest.approx([-1.5902428721, -0.6627677745, 0.1437164356], rel=1e-8)
    assert r["p"][[0, 6, 22]] == pytest.approx([0.1368880684, 0.5195845657, 0.8880289056], rel=1e-8)


def test_bh_and_half_sample_mode():
    from scipy.stats import false_discovery_control

    p = np.random.default_rng(0).uniform(size=200) ** 3
    assert bh(p) == pytest.approx(false_discovery_control(p))
    # 30 % of values shifted one way: the mode stays with the majority, the median does not.
    x = np.concatenate([np.random.default_rng(1).normal(0, 0.1, 700), np.random.default_rng(2).normal(1.0, 0.1, 300)])
    assert abs(half_sample_mode(x)) < 0.03 and np.median(x) > 0.05


def test_site_specific_changes_found_under_global_scaling():
    profiles, truth = simulate("global_with_site", SimSettings(), seed=17)
    a = analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=10))
    rel = per_region_tests(a).tests["log2_signal_relative"]
    changed = truth["regions"]["amp_rel"][rel["region_index"]] == 1
    called = rel["call"] == "Higher"
    assert called[changed].mean() >= 0.9  # power
    false = np.isin(rel["call"], ("Higher", "Lower")) & ~changed
    assert false.sum() <= 0.1 * max(called.sum(), 1)  # FDR-like bound


def test_uniform_efficiency_change_is_not_site_specific():
    # Peak-only scaling changes window AUC less in background-dominated regions; the relative metric
    # is built on background-corrected signal so it must stay quiet.
    a = analyse(simulate("amplitude", SimSettings(), seed=19)[0], Config(reference="WT", alternative="KO", bootstrap=10))
    tests = per_region_tests(a).tests
    assert np.isin(tests["log2_signal_relative"]["call"], ("Higher", "Lower")).sum() <= 2
    assert (tests["log2_signal"]["call"] == "Higher").mean() > 0.5


def test_per_region_shift_finds_shifts_in_both_directions():
    profiles, truth = simulate("opposite_shift", SimSettings(), seed=23)
    a = analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=10))
    res = per_region_tests(a)
    sh = res.tests["shift"]
    direction = truth["regions"]["shift"][sh["region_index"]]
    assert (sh["call"][direction > 0] == "Higher").mean() > 0.6
    assert (sh["call"][direction < 0] == "Lower").mean() > 0.6
    assert np.median(sh["estimate"][direction > 0]) > 50  # template fit: little attenuation
    assert (res.tests["shape_kernel"]["call"] == "Changed").mean() > 0.6
    robust = {r["metric"]: r for r in res.robust_global}
    assert abs(robust["median_region_centroid"]["estimate"]) < 30  # no net shift


def test_kernel_shape_test_two_vs_two_and_amplitude_only():
    """2 vs 2 works (assignments pooled over regions); amplitude-only changes are not shape calls."""
    for scenario, n in (("shift", 2), ("amplitude", 3), ("site_amplitude", 2)):
        profiles, _ = simulate(scenario, SimSettings(n_ref=n, n_alt=n), seed=41)
        res = per_region_tests(analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=10)))
        changed = (res.tests["shape_kernel"]["call"] == "Changed").mean()
        if scenario == "shift":
            assert changed > 0.5
        else:
            assert changed < 0.01


def test_shape_classes_are_condition_blind_and_alignment_recovers_anchor_error():
    from genomecanvas.profileseg.data import Profiles
    from genomecanvas.profileseg.shapeclass import align_profiles, pooled_core, shape_classes

    prof, _ = simulate("null", SimSettings(n_regions=400), seed=5)
    rng = np.random.default_rng(2)
    jit = rng.integers(-5, 6, 400)
    v = prof.values.copy()
    for r, q in enumerate(jit):
        v[:, r] = np.roll(v[:, r], q, axis=-1)
    prof = Profiles(v, prof.edges, prof.samples, prof.regions, prof.oriented)
    P, ok = pooled_core(prof.values, prof.edges, 500, 1000)
    sc = shape_classes(P[ok], 8)
    assert np.corrcoef(sc.shift_bins, -jit[ok])[0, 1] > 0.8
    # relabelling the samples does not change the pooled profiles, hence classes and shifts
    swapped = Profiles(v[::-1], prof.edges, prof.samples, prof.regions, prof.oriented)
    P2, _ = pooled_core(swapped.values, swapped.edges, 500, 1000)
    assert np.allclose(P, P2)
    aligned, shift_bp, _ = align_profiles(prof, 500, 1000, 100.0)
    assert np.all(np.abs(shift_bp) <= 100.0) and aligned.values.shape == prof.values.shape


def test_stratified_pooled_permutation_p_is_uniform_under_exchangeability():
    from genomecanvas.profileseg.regional import pooled_permutation_p

    rng = np.random.default_rng(0)
    # two strata with different null scales; rows are exchangeable within each region
    F = np.concatenate([rng.exponential(1.0, (10, 500)), rng.exponential(5.0, (10, 500))], axis=1)
    strata = np.r_[np.zeros(500, int), np.ones(500, int)]
    p = pooled_permutation_p(F, 0, strata)
    assert abs(np.mean(p < 0.05) - 0.05) < 0.02
    assert abs(np.mean(p[strata == 1] < 0.05) - 0.05) < 0.03  # without strata this stratum would be liberal
    p_global = pooled_permutation_p(F, 0, None)
    assert np.mean(p_global[strata == 1] < 0.05) > 0.075


def test_shape_class_summaries_are_descriptive_and_optional():
    profiles, _ = simulate("shift", SimSettings(n_regions=300), seed=8)
    off = per_region_tests(analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=10)))
    on = per_region_tests(analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=10, shape_classes=True)))
    assert off.classes == []
    assert on.classes and all("median_abs_shift_called_bp" in r for r in on.classes)
    # classes never change the tests
    assert np.array_equal(off.tests["shape_kernel"]["p"], on.tests["shape_kernel"]["p"])
    assert np.array_equal(off.tests["shift"]["p"], on.tests["shift"]["p"])


def test_kernel_feature_map_matches_explicit_kernel_quadratic_forms():
    """kernel_shape_test uses φ = pL with K = LLᵀ; the statistics must equal the explicit xᵀKx forms."""
    import itertools

    from genomecanvas.profileseg.regional import kernel_shape_test, pooled_permutation_p, support_mask

    rng = np.random.default_rng(4)
    y = rng.gamma(2.0, 1.0, (6, 200, 40))
    c = np.arange(40) * 20.0
    ref, alt = np.arange(3), np.arange(3, 6)
    p = y / (y * support_mask(y, 0.2)[None]).sum(-1, keepdims=True)  # same normalisation as the test
    K = np.exp(-0.5 * ((c[:, None] - c[None]) / 60.0) ** 2)
    dev = p - p.mean(0)
    s0 = np.quantile(np.einsum("srb,bc,src->r", dev, K, dev) / 5, 0.5)
    splits = [s for s in itertools.combinations(range(6), 3) if 0 in s]
    F = []
    for s in splits:
        a = np.array(s); b = np.setdiff1d(np.arange(6), a)
        d = p[a].mean(0) - p[b].mean(0)
        res = np.concatenate([p[a] - p[a].mean(0), p[b] - p[b].mean(0)])
        F.append(np.einsum("rb,bc,rc->r", d, K, d) / (np.einsum("srb,bc,src->r", res, K, res) / 4 + s0))
    expected = pooled_permutation_p(np.array(F), 0)
    got = kernel_shape_test(y, np.ones((200, 40), bool), ref, alt, 60.0, 0.5, c)
    assert np.allclose(got, expected)


def test_batch_confounded_design_is_not_tested_and_within_batch_assignments():
    from genomecanvas.profileseg.regional import label_assignments

    # 2 vs 2 with one replicate of each condition per batch: two distinct assignments within batches
    splits, obs = label_assignments(2, 4, np.array(["b1", "b2", "b1", "b2"]))
    assert len(splits) == 2 and splits[obs] == (0, 1)
    # condition = batch: only the observed assignment
    assert len(label_assignments(2, 4, np.array(["b1", "b1", "b2", "b2"]))[0]) == 1
    profiles, _ = simulate("shift", SimSettings(n_regions=200, n_ref=2, n_alt=2), seed=9)
    for k, smp in enumerate(profiles.samples):
        smp["batch"] = "b1" if k < 2 else "b2"
    res = per_region_tests(analyse(profiles, Config(reference="WT", alternative="KO", bootstrap=10)))
    assert np.all(res.tests["shape_kernel"]["call"] == "not_tested")
    assert any("condition and batch coincide" in n for n in res.notes)


def test_kernel_test_does_not_turn_amplitude_changes_into_shape_calls():
    """Signed normalisation: with 4 vs 4 and fine bins, clipping made amplitude-only changes look like broadening."""
    st = SimSettings(n_ref=4, n_alt=4, flank=3000, binsize=10, n_regions=250)
    for scenario in ("amplitude", "site_amplitude"):
        profiles, _ = simulate(scenario, st, seed=31)
        cfg = Config(reference="WT", alternative="KO", flank=3000, binsize=10, bootstrap=10)
        res = per_region_tests(analyse(profiles, cfg, bootstrap=0))
        assert (res.tests["shape_kernel"]["call"] == "Changed").mean() < 0.005
