#!/usr/bin/python
"""Functions related to distance computations between distributions."""

from typing import Any
import torch
import numpy as np
from torch_cluster import knn
from scipy.stats import wasserstein_distance

# LOCAL IMPORTS
from .common_asserts import check_type, check_num_value


########################################


def _get_distance_closest(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Computes the smallest distance between two structures as coordinate arrays."""
    src: torch.Tensor = torch.from_numpy(source) # type: ignore
    tgt: torch.Tensor = torch.from_numpy(target) # type: ignore

    idx_src, idx_tgt = knn(tgt, src, 1)

    closest_distance = (src[idx_src] - tgt[idx_tgt]).norm(dim=1)

    return closest_distance.numpy()


########################################


def recall(source: np.ndarray, target: np.ndarray, threshold: float) -> float:
    """Get Recall metric value."""
    check_type(source, "source", (np.ndarray,))
    check_type(target, "target", (np.ndarray,))
    check_type(threshold, "threshold", (float,))
    check_num_value(threshold, "threshold", ">", 0)

    distance = _get_distance_closest(source=target, target=source)

    mask = distance < threshold

    return mask.astype(np.float32).mean().item()


########################################


def precision(source: np.ndarray, target: np.ndarray, threshold: float) -> float:
    """Get Precision metric value."""
    check_type(source, "source", (np.ndarray,))
    check_type(target, "target", (np.ndarray,))
    check_type(threshold, "threshold", (float,))
    check_num_value(threshold, "threshold", ">", 0)

    distance = _get_distance_closest(source=source, target=target)

    mask = distance < threshold

    return mask.astype(np.float32).mean().item()


########################################


def _sqrt_tr(x: np.ndarray) -> float:
    """Trace of the square root eigenvalues of a matrix."""
    return np.sum(np.sqrt(np.linalg.eigvals(x).real))


########################################


def frechet_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Frechet Distance between two matrices x and y."""
    check_type(x, "x", (np.ndarray,))
    check_type(y, "y", (np.ndarray,))

    mu_x = np.mean(x, axis=0)
    mu_y = np.mean(y, axis=0)

    sigma_x = np.cov(x.T)
    sigma_y = np.cov(y.T)

    fid = (
        np.linalg.norm(mu_x - mu_y) ** 2
        + sigma_x.trace()
        + sigma_y.trace()
        - 2 * _sqrt_tr(sigma_x @ sigma_y)
    )

    return fid.item()


########################################


def get_emd(*args, **kwargs) -> Any:
    """Simple wrapper for scipy.stats.wasserstein_distance."""
    return wasserstein_distance(*args, **kwargs)


########################################
