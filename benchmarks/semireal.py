"""Semi-real data: real replicate noise from GSE112560 bigWigs around a known truth.

Seven real samples (ZNF384 WT ×3, TAF15 ×2, TCF3 ×2) over the same consensus regions. After
per-sample scaling, Helmert contrasts within each condition give four independent real noise
vectors (WT 2, TAF15 1, TCF3 1): spatially correlated, heteroscedastic, region-specific. A
pseudo-experiment takes the WT mean profile as the base, applies a known change to the base of
the second group, and adds one noise vector per pseudo-sample (rescaled for the base level by
√(base / condition mean), a Poisson-like variance law). 2 vs 2 uses four distinct vectors at the
same regions. Larger designs reuse vectors, taken from donor regions of similar strength so that
pseudo-samples stay independent (region-specific noise structure is then only approximate).

An earlier version used plain residuals from the condition means; with two samples those are
exactly r and −r, so pseudo-samples were dependent. Results from that version are superseded.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from genomecanvas.profileseg.data import Profiles, extract, load_config, read_manifest, read_regions


def build_cache(config_path: Path, out: Path) -> Path:
    config = load_config(config_path)
    samples, _ = read_manifest(config.path(config.manifest))
    regions = read_regions(config.path(config.regions))
    prof = extract(config, samples, regions)
    np.savez_compressed(out, values=prof.values, edges=prof.edges,
                        conditions=np.array([s["condition"] for s in samples]))
    return out


class SemiReal:
    def __init__(self, cache: Path, min_mass_quantile: float = 0.2):
        d = np.load(cache)
        v, self.edges, cond = d["values"], d["edges"], d["conditions"]
        keep = np.all(np.isfinite(v), axis=(0, 2))
        v = v[:, keep]
        v = v / np.median(v.sum(-1), axis=1)[:, None, None] * np.median(v.sum(-1))  # per-sample scaling
        self.base = v[cond == "WT"].mean(0)
        # Independent noise vectors: Helmert contrasts within each condition (n samples give n − 1
        # orthonormal contrasts, each with the replicate variance). Plain residuals from the mean are
        # dependent (with two samples they are exactly r and −r), which breaks exchangeability.
        res, level = [], []
        for c in np.unique(cond):
            x = v[cond == c]
            m = x.mean(0)
            for j in range(1, len(x)):
                h = (x[:j].sum(0) - j * x[j]) / np.sqrt(j * (j + 1))
                res.append(h)
                level.append(m)
        self.res, self.level = np.array(res), np.array(level)
        mass = self.base.sum(-1)
        self.pool = np.flatnonzero(mass >= np.quantile(mass, min_mass_quantile))
        # mass deciles for donor regions (a noise vector reused for a second sample is taken from
        # another region of similar strength, so that pseudo-samples stay independent)
        self.decile = np.clip((np.argsort(np.argsort(mass)) * 10) // len(mass), 0, 9)

    def experiment(self, change: str, n_regions: int, fraction: float, seed: int, n_per_group: int = 3,
                   shift_bp: float = 100.0, broaden: float = 1.5, fold: float = 2.0) -> tuple[Profiles, np.ndarray]:
        rng = np.random.default_rng(seed)
        idx = rng.choice(self.pool, n_regions, replace=False)
        base = self.base[idx]
        c = 0.5 * (self.edges[:-1] + self.edges[1:])
        changed = np.zeros(n_regions, bool)
        changed[: int(round(fraction * n_regions))] = change != "null"
        alt = base.copy()
        for r in np.flatnonzero(changed):
            if change == "shift":
                alt[r] = np.interp(c - shift_bp, c, base[r], left=base[r][0], right=base[r][-1])
            elif change == "broaden":
                alt[r] = np.interp(c / broaden, c, base[r]) * 1.0
                alt[r] *= base[r].sum() / max(alt[r].sum(), 1e-9)
            elif change == "amplitude":  # region-specific amplitude, no shape change
                bg = np.r_[base[r][:25], base[r][-25:]].mean()
                alt[r] = bg + (base[r] - bg) * fold
        if change == "global_scale":
            alt = alt * 0.6
            changed[:] = False
        n = n_per_group
        K = len(self.res)
        pick = np.concatenate([rng.permutation(K) for _ in range(-(-2 * n // K))])[: 2 * n]
        used = {}
        values = np.empty((2 * n, n_regions, base.shape[1]))
        for s, k in enumerate(pick):
            mean = base if s < n else alt
            if k in used:  # reuse of a noise vector: draw it from matched donor regions instead
                donor = np.array([rng.choice(self.pool[self.decile[self.pool] == self.decile[r]]) for r in idx])
            else:
                donor = idx
            used[k] = True
            lvl = self.level[k][donor]
            scale = np.sqrt((mean.mean(-1, keepdims=True) + 1e-3) / (lvl.mean(-1, keepdims=True) + 1e-3)) if k in used and donor is not idx \
                else np.sqrt((mean + 1e-3) / (lvl + 1e-3))
            values[s] = np.clip(mean + self.res[k][donor] * np.clip(scale, 0, 5), 0, None)
        samples = [dict(sample_id=f"{'WT' if s < n else 'KO'}_rep{s % n + 1}", condition="WT" if s < n else "KO",
                        biological_replicate=str(s % n + 1)) for s in range(2 * n)]
        regions = [dict(chrom="chrSemi", start=int(10_000 * (i + 1)), end=int(10_000 * (i + 1) + 1), strand="+",
                        name=f"r{i}") for i in range(n_regions)]
        truth_changed = changed if change in ("shift", "broaden") else np.zeros(n_regions, bool)
        return Profiles(values, self.edges.copy(), samples, regions, True), truth_changed


if __name__ == "__main__":
    import sys
    print(build_cache(Path(sys.argv[1]), Path(sys.argv[2])))
