"""Replicate-aware profile simulation with known truth (parametric or real-data templates)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .data import Profiles

# Scenario → which estimand truly changes (used by the benchmark to score calls).
SCENARIOS = {
    "null": dict(truth={}),
    "amplitude": dict(
        truth={"log2_auc": +1, "log2_centre_auc": +1}
    ),  # peak ×1.6 over constant background
    "global_scale": dict(truth={"log2_auc": +1, "log2_centre_auc": +1}),  # peak and background ×1.6
    "shift": dict(truth={"centroid": +1, "asymmetry": +1}),  # +100 bp translation
    "broaden": dict(truth={"sd_width": +1, "iqr_width": +1}),  # σ ×1.5 at conserved area
    "right_flank": dict(truth={"asymmetry": +1, "centroid": +1, "log2_auc": +1, "sd_width": +1}),
    "left_flank": dict(truth={"asymmetry": -1, "centroid": -1, "log2_auc": +1, "sd_width": +1}),
    "sparse_shift": dict(truth={"centroid": +1, "asymmetry": +1}),  # 25 % of regions shifted
    "fragment_length": dict(truth={"sd_width": +1, "iqr_width": +1}),  # technical: wider fragments
    # Region-heterogeneous scenarios (aggregate truth comes from the noise-free profiles):
    "opposite_shift": dict(truth={}),  # half the regions +100 bp, half −100 bp
    "site_amplitude": dict(truth={}),  # 20 % of regions gain ×2.5 (site-specific)
    "global_with_site": dict(truth={}),  # everything ×0.6, and 20 % of regions ×2.5 on top
}


@dataclass
class SimSettings:
    n_ref: int = 3
    n_alt: int = 3
    n_regions: int = 300
    flank: int = 2000
    binsize: int = 20
    depth: float = 40.0  # expected reads per bin at the peak apex for a typical region
    background: float = 0.25  # background relative to peak apex
    dispersion: float = 0.15  # negative-binomial dispersion (variance = μ + φμ²)
    replicate_sd: float = 0.15  # log-scale replicate amplitude noise per region
    library_sd: float = 0.10  # log-scale library-size noise per sample
    ar_rho: float = 0.6  # autocorrelation of multiplicative bin noise
    ar_sd: float = 0.08
    shift_bp: float = 100.0
    broaden: float = 1.5
    amplitude: float = 1.6
    sparse_fraction: float = 0.25
    site_fraction: float = 0.2  # share of regions with a site-specific gain
    site_fold: float = 2.5
    fragment_extra_bp: float = 120.0
    oriented: bool = True


def _gauss(c, mu, sigma):
    return np.exp(-0.5 * ((c - mu) / sigma) ** 2)


def _smooth(profile: np.ndarray, sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0:
        return profile
    half = int(np.ceil(4 * sigma_bins))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) / sigma_bins) ** 2)
    k /= k.sum()
    return np.apply_along_axis(lambda r: np.convolve(r, k, mode="same"), -1, profile)


def _box(profile: np.ndarray, width_bins: int) -> np.ndarray:
    if width_bins <= 1:
        return profile
    k = np.ones(width_bins) / width_bins
    return np.apply_along_axis(lambda r: np.convolve(r, k, mode="same"), -1, profile)


def base_shapes(s: SimSettings, rng, templates: np.ndarray | None = None):
    """Per-region peak component (apex ≈ 1) and background level, before any scenario."""
    edges = np.arange(-s.flank, s.flank + s.binsize, s.binsize, dtype=float)
    c = (edges[:-1] + edges[1:]) / 2
    if templates is not None:
        pick = rng.choice(len(templates), s.n_regions, replace=len(templates) < s.n_regions)
        t = np.asarray(templates, dtype=float)[pick]
        if t.shape[1] != len(c):
            raise ValueError("template bins do not match the simulation grid")
        n = max(1, len(c) // 8)
        bg = np.median(np.concatenate([t[:, :n], t[:, -n:]], axis=1), axis=1)
        peak = np.clip(t - bg[:, None], 0, None)
        apex = peak.max(axis=1, keepdims=True)
        peak = np.where(apex > 0, peak / np.where(apex > 0, apex, 1), 0)
        bg_rel = np.clip(bg / np.maximum(apex[:, 0], 1e-12), 0.05, 2.0)
    else:
        mu = rng.normal(0, 30, s.n_regions)
        sigma = rng.lognormal(np.log(150), 0.25, s.n_regions)
        peak = _gauss(c[None, :], mu[:, None], sigma[:, None])
        bg_rel = np.full(s.n_regions, s.background)
    strength = rng.lognormal(0, 0.5, s.n_regions)  # between-region heterogeneity
    return edges, c, peak, bg_rel, strength


def apply_scenario(
    name: str, peak: np.ndarray, bg_rel: np.ndarray, c: np.ndarray, s: SimSettings, rng
):
    """Transform the alternative condition's peak/background; returns per-region truth too.

    Region truth (sign per region): ``amp_abs`` absolute amplitude change, ``amp_rel`` change
    relative to the bulk of regions (site-specific), ``shift`` direction of translation.
    """
    peak, bg = peak.copy(), bg_rel.copy()
    step = s.binsize
    n = len(peak)
    truth = dict(amp_abs=np.zeros(n, int), amp_rel=np.zeros(n, int), shift=np.zeros(n, int))
    if name in ("null",):
        pass
    elif name == "amplitude":
        peak = peak * s.amplitude
        truth["amp_abs"][:] = 1
    elif name == "global_scale":
        peak, bg = peak * s.amplitude, bg * s.amplitude
        truth["amp_abs"][:] = 1
    elif name in ("site_amplitude", "global_with_site"):
        rows = rng.choice(n, int(s.site_fraction * n), replace=False)
        if name == "global_with_site":
            peak, bg = peak * 0.6, bg * 0.6
            truth["amp_abs"][:] = -1
        peak[rows] *= s.site_fold
        truth["amp_abs"][rows] = 1
        truth["amp_rel"][rows] = 1
    elif name == "opposite_shift":
        for r in range(n):
            sign = 1 if r % 2 else -1
            peak[r] = np.interp(c - sign * s.shift_bp, c, peak[r], left=0, right=0)
            truth["shift"][r] = sign
    elif name in ("shift", "sparse_shift"):
        rows = (
            np.arange(len(peak))
            if name == "shift"
            else rng.choice(len(peak), int(s.sparse_fraction * len(peak)), replace=False)
        )
        for r in rows:
            peak[r] = np.interp(c - s.shift_bp, c, peak[r], left=0, right=0)
            truth["shift"][r] = 1
    elif name == "broaden":
        # Convolution conserves area; extra σ so that σ_new ≈ broaden × σ_old for a ~150 bp peak.
        extra = 150 * np.sqrt(s.broaden**2 - 1) / step
        peak = _smooth(peak, extra)
    elif name in ("right_flank", "left_flank"):
        sign = 1 if name == "right_flank" else -1
        peak = peak + 0.5 * _gauss(c[None, :], sign * 350, 120) * peak.max(axis=1, keepdims=True)
    elif name == "fragment_length":
        half = int(round(s.fragment_extra_bp / step / 2))
        peak = _box(peak, 2 * half + 1)  # odd width keeps the kernel centred (no artificial shift)
    else:
        raise ValueError(f"unknown scenario {name}; choose from {sorted(SCENARIOS)}")
    return peak, bg, truth


def _ar_noise(rng, shape, rho, sd):
    e = rng.normal(0, sd, shape)
    out = np.empty(shape)
    out[..., 0] = e[..., 0]
    for k in range(1, shape[-1]):
        out[..., k] = rho * out[..., k - 1] + np.sqrt(1 - rho**2) * e[..., k]
    return np.exp(out - sd**2 / 2)


def simulate(
    scenario: str,
    settings: SimSettings | None = None,
    seed: int = 0,
    templates: np.ndarray | None = None,
    noiseless: bool = False,
) -> tuple[Profiles, dict]:
    """Simulate a scenario. ``noiseless`` returns the expected profiles (for computing true effects)."""
    s = settings or SimSettings()
    rng = np.random.default_rng(seed)
    edges, c, peak, bg_rel, strength = base_shapes(s, rng, templates)
    peak_alt, bg_alt, region_truth = apply_scenario(scenario, peak, bg_rel, c, s, rng)
    n = s.n_ref + s.n_alt
    values = np.empty((n, s.n_regions, len(c)))
    samples = []
    libraries = []
    for i in range(n):
        alt = i >= s.n_ref
        pk, bg = (peak_alt, bg_alt) if alt else (peak, bg_rel)
        # Mean-one log-normal factors, so the noise-free profile is the expectation of the noisy one.
        library = rng.lognormal(-(s.library_sd**2) / 2, s.library_sd)
        libraries.append(float(library))
        rep = rng.lognormal(-(s.replicate_sd**2) / 2, s.replicate_sd, s.n_regions)[:, None]
        if noiseless:
            values[i] = s.depth * strength[:, None] * (pk + bg[:, None])
        else:
            mean = s.depth * library * strength[:, None] * (rep * pk + bg[:, None])
            mean = mean * _ar_noise(rng, mean.shape, s.ar_rho, s.ar_sd)
            r = 1.0 / s.dispersion
            counts = rng.negative_binomial(r, r / (r + mean))
            values[i] = counts / library  # library-size normalised coverage (floating point)
        cond = "KO" if alt else "WT"
        samples.append(
            dict(
                sample_id=f"{cond}_rep{(i - s.n_ref if alt else i) + 1}",
                condition=cond,
                biological_replicate=str((i - s.n_ref if alt else i) + 1),
                signal_type="simulated",
            )
        )
    regions = [
        dict(
            chrom="chrSim",
            start=10_000 + 10_000 * k,
            end=10_001 + 10_000 * k,
            name=f"sim{k}",
            strand="+" if s.oriented else ".",
        )
        for k in range(s.n_regions)
    ]
    truth = dict(
        scenario=scenario,
        changes=SCENARIOS[scenario]["truth"],
        shift_bp=s.shift_bp
        if scenario == "shift"
        else (s.shift_bp * s.sparse_fraction if scenario == "sparse_shift" else 0.0),
        amplitude=s.amplitude if scenario in ("amplitude", "global_scale") else 1.0,
        regions=region_truth,
        library=libraries,  # values × library = integer counts (for count-model comparators)
    )
    return Profiles(
        values,
        edges,
        samples,
        regions,
        s.oriented,
        provenance=dict(source="simulation", seed=seed, scenario=scenario),
    ), truth


def write_example(
    directory, settings: SimSettings | None = None, seed: int = 7201, scenario: str = "shift"
) -> dict:
    """Synthetic bigWigs + BED + manifest + config, so the example exercises real extraction."""
    from pathlib import Path

    import pyBigWig

    s = settings or SimSettings(n_regions=200)
    directory = Path(directory)
    (directory / "bigwig").mkdir(parents=True, exist_ok=True)
    (directory / "metadata").mkdir(exist_ok=True)
    profiles, truth = simulate(scenario, s, seed)
    chrom_len = 10_000 * (s.n_regions + 2)
    n_bins = profiles.values.shape[2]
    rows = []
    for i, sample in enumerate(profiles.samples):
        path = directory / "bigwig" / f"{sample['sample_id']}.bw"
        bw = pyBigWig.open(str(path), "w")
        bw.addHeader([("chrSim", chrom_len)])
        for k, region in enumerate(profiles.regions):
            anchor = region["start"]
            starts = anchor - s.flank + s.binsize * np.arange(n_bins)
            bw.addEntries(
                ["chrSim"] * n_bins,
                starts.tolist(),
                ends=(starts + s.binsize).tolist(),
                values=profiles.values[i, k].astype(float).tolist(),
            )
        bw.close()
        rows.append(
            {
                **sample,
                "assay": "ChIP-seq",
                "genome_build": "simulated",
                "signal_file": f"bigwig/{path.name}",
                "normalization": "library size",
                "region_set": "sim",
                "batch": "1",
                "pairing_id": "",
                "bigwig_params": f"binSize={s.binsize}",
                "fragment_length": "200",
            }
        )
    with open(directory / "regions.bed", "w") as handle:
        for region in profiles.regions:
            handle.write(
                f"{region['chrom']}\t{region['start']}\t{region['end']}\t{region['name']}\t0\t{region['strand']}\n"
            )
    from .data import MANIFEST_COLUMNS

    with open(directory / "metadata" / "sample_manifest.tsv", "w") as handle:
        handle.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in MANIFEST_COLUMNS) + "\n")
    (directory / "config.yaml").write_text(
        "manifest: metadata/sample_manifest.tsv\nregions: regions.bed\n"
        "contrast: {reference: WT, alternative: KO}\nanchor: start\norientation: auto\n"
        f"flank: {s.flank}\nbinsize: {s.binsize}\nbackground: {{width: 500}}\n"
        "shape: {core_flank: 1000, negative_policy: signed, min_mass_quantile: 0.1}\n"
        "inference: {bootstrap: 1000, seed: 7201}\noutput: results\n",
        encoding="utf-8",
    )
    return truth
