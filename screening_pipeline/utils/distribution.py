"""Functions related to distance computations between distributions."""

from typing import Any
import torch
import numpy as np
from torch_cluster import knn
from scipy.stats import wasserstein_distance


def get_distance_closest(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Computes the smallest distance between two structures as coordinate arrays."""
    source = torch.from_numpy(source)
    target = torch.from_numpy(target)

    idx_src, idx_tgt = knn(target, source, 1)

    closest_distance = (source[idx_src] - target[idx_tgt]).norm(dim=1)

    return closest_distance.numpy()


def recall(source: np.ndarray, target: np.ndarray, threshold: float) -> float:
    """Get Recall metric value."""
    distance = get_distance_closest(source=target, target=source)

    mask = distance < threshold

    return mask.astype(np.float32).mean().item()


def precision(source: np.ndarray, target: np.ndarray, threshold: float) -> float:
    """Get Precision metric value."""
    distance = get_distance_closest(source=source, target=target)

    mask = distance < threshold

    return mask.astype(np.float32).mean().item()


def sqrt_tr(x: np.ndarray) -> float:
    """Trace of the square root eigenvalues of a matrix."""
    return np.sum(np.sqrt(np.linalg.eigvals(x).real))


def frechet_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Frechet Distance between two matrices x and y."""
    mu_x = np.mean(x, axis=0)
    mu_y = np.mean(y, axis=0)

    sigma_x = np.cov(x.T)
    sigma_y = np.cov(y.T)

    fid = (
        np.linalg.norm(mu_x - mu_y) ** 2
        + sigma_x.trace()
        + sigma_y.trace()
        - 2 * sqrt_tr(sigma_x @ sigma_y)
    )

    return fid.item()


def emd_wrapper(*args, **kwargs) -> Any:
    """Simple wrapper for scipy.stats.wasserstein_distance."""
    return wasserstein_distance(*args, **kwargs)