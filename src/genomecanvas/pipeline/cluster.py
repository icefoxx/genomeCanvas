"""Region clustering on an enrichment matrix."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

from .matrix import EnrichmentMatrix, compute_region_matrix
from .metadata import SampleMeta

ClusterMethod = Literal[
    "none", "kmeans", "hclust", "kmeansExt", "hclustExt", "miniBatchKMeansExt",
    "leidenExt", "leidenKmeansExt",
]


def cluster_labels(
    matrix: EnrichmentMatrix,
    *,
    method: ClusterMethod,
    n_clusters: int = 5,
    sample_indices: list[int] | None = None,
    dim: str = "raw",
    seed: int | None = 7201,
    leiden_resolution: float = 1.0,
) -> np.ndarray:  # type: ignore[type-arg]
    """Return 1-based cluster ids for each region."""
    if method in {"none"}:
        return np.ones(matrix.values.shape[0], dtype=int)
    features = _feature_matrix(matrix, sample_indices)
    space = _embed(features, dim, seed)
    if method in {"kmeans", "kmeansExt"}:
        return _kmeans(space, n_clusters, seed)
    if method == "miniBatchKMeansExt":
        return _kmeans(space, n_clusters, seed, mini_batch=True)
    if method in {"hclust", "hclustExt"}:
        return _hclust(space, n_clusters)
    if method == "leidenExt":
        return _leiden(space, leiden_resolution, seed)
    if method == "leidenKmeansExt":
        leiden = _leiden(space, leiden_resolution, seed)
        unique = np.unique(leiden)
        if n_clusters >= unique.size:
            return leiden
        centroids = np.vstack([space[leiden == label].mean(axis=0) for label in unique])
        merged = _kmeans(centroids, n_clusters, seed)
        mapping = {int(label): int(merged[index]) for index, label in enumerate(unique)}
        return np.array([mapping[int(label)] for label in leiden], dtype=int)
    raise ValueError(f"unsupported clustering method {method}")


def write_cluster_beds(
    matrix: EnrichmentMatrix, labels: np.ndarray, prefix: str | Path  # type: ignore[type-arg]
) -> list[tuple[str, str]]:
    """Write one BED per cluster and return (path, label) pairs."""
    destination = Path(prefix)
    destination.parent.mkdir(parents=True, exist_ok=True)
    beds: list[tuple[str, str]] = []
    for cluster in np.unique(labels):
        path = Path(f"{destination}.cluster{int(cluster)}.bed")
        rows = [
            f"{coord.chrom}\t{coord.start}\t{coord.end}"
            for coord, label in zip(matrix.coords, labels, strict=True)
            if label == cluster
        ]
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        beds.append((str(path), f"C{int(cluster)}"))
    return beds


def recompute_cluster_matrix(
    samples: list[SampleMeta],
    beds: list[tuple[str, str]],
    *,
    flank: int,
    binsize: int,
    ref_point: str,
    body_length: int,
    engine: str,
) -> EnrichmentMatrix:
    """Rebuild an enrichment matrix using cluster BEDs as groups."""
    return compute_region_matrix(
        samples,
        beds,
        flank=flank,
        binsize=binsize,
        ref_point=ref_point,  # type: ignore[arg-type]
        body_length=body_length,
        engine=engine,  # type: ignore[arg-type]
    )


def _feature_matrix(
    matrix: EnrichmentMatrix, sample_indices: list[int] | None
) -> np.ndarray:  # type: ignore[type-arg]
    indices = sample_indices or list(range(len(matrix.sample_labels)))
    blocks = [matrix.sample_block(index) for index in indices]
    joined = np.hstack(blocks) if blocks else matrix.values
    return np.nan_to_num(joined, nan=0.0)


def _embed(features: np.ndarray, dim: str, seed: int | None) -> np.ndarray:  # type: ignore[type-arg]
    if dim == "raw":
        return features
    centered = features - features.mean(axis=0, keepdims=True)
    _u, _s, vt = np.linalg.svd(centered, full_matrices=False)
    scores = centered @ vt[: min(50, vt.shape[0])].T
    if dim == "pca":
        return scores
    if dim == "tsne":
        try:
            from sklearn.manifold import TSNE
        except ImportError as error:
            raise ImportError("t-SNE clustering requires scikit-learn") from error
        perplexity = min(30, max(5, features.shape[0] // 4))
        embedded = TSNE(
            n_components=2, init="pca", random_state=seed, perplexity=perplexity
        ).fit_transform(scores[:, : min(50, scores.shape[1])])
        return embedded
    if dim == "umap":
        try:
            from umap import UMAP
        except ImportError as error:
            raise ImportError("UMAP clustering requires umap-learn") from error
        reduced = scores[:, : min(50, scores.shape[1])]
        return UMAP(n_components=2, random_state=seed).fit_transform(reduced)
    raise ValueError("clustering dim must be raw, pca, tsne, or umap")


def _kmeans(
    space: np.ndarray, n_clusters: int, seed: int | None, *, mini_batch: bool = False  # type: ignore[type-arg]
) -> np.ndarray:  # type: ignore[type-arg]
    try:
        from sklearn.cluster import KMeans, MiniBatchKMeans
    except ImportError:
        return _numpy_kmeans(space, n_clusters, seed)
    if mini_batch:
        model = MiniBatchKMeans(
            n_clusters=n_clusters, random_state=seed, n_init=10, batch_size=4096
        )
    else:
        model = KMeans(n_clusters=n_clusters, random_state=seed, n_init=25)
    return model.fit_predict(space) + 1


def _numpy_kmeans(space: np.ndarray, n_clusters: int, seed: int | None) -> np.ndarray:  # type: ignore[type-arg]
    rng = np.random.default_rng(seed)
    centers = space[rng.choice(space.shape[0], n_clusters, replace=False)]
    labels = np.zeros(space.shape[0], dtype=int)
    for _ in range(25):
        distances = ((space[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        labels = distances.argmin(axis=1)
        for cluster in range(n_clusters):
            members = space[labels == cluster]
            if members.size:
                centers[cluster] = members.mean(axis=0)
    return labels + 1


def _hclust(space: np.ndarray, n_clusters: int) -> np.ndarray:  # type: ignore[type-arg]
    try:
        from scipy.cluster.hierarchy import fcluster, linkage
    except ImportError as error:
        raise ImportError("hierarchical clustering requires scipy") from error
    tree = linkage(space, method="ward")
    return fcluster(tree, t=n_clusters, criterion="maxclust").astype(int)


def _leiden(space: np.ndarray, resolution: float, seed: int | None) -> np.ndarray:  # type: ignore[type-arg]
    try:
        import igraph as ig
        import leidenalg as la
        from sklearn.neighbors import NearestNeighbors
    except ImportError as error:
        raise ImportError(
            "leiden clustering requires igraph, leidenalg, and scikit-learn"
        ) from error
    n_neighbors = max(5, min(15, space.shape[0] - 1))
    neighbors = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean").fit(space)
    distances, indices = neighbors.kneighbors(space)
    edges, weights = [], []
    for i in range(space.shape[0]):
        for dist, j in zip(distances[i, 1:], indices[i, 1:], strict=True):
            if i < int(j):
                edges.append((i, int(j)))
                weights.append(1.0 / (1.0 + float(dist)))
    graph = ig.Graph(n=space.shape[0], edges=edges, directed=False)
    graph.es["weight"] = weights
    partition = la.find_partition(
        graph,
        la.RBConfigurationVertexPartition,
        weights="weight",
        seed=seed,
        resolution_parameter=float(resolution),
    )
    return np.array(partition.membership, dtype=int) + 1
