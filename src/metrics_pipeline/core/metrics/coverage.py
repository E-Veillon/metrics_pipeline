"""Compute Coverage (Precision) and Coverage (Recall) metrics."""

import typing as tp
import functools as ft

import numpy as np

import torch
from torch_cluster import knn

from pymatgen.core import Structure

from .metric_base import StructureDistribution, Metric, MetricsData
from core.utils import GenMatStructure


class Coverage(Metric):
    """
    Compute Coverage (Precision) and Coverage (Recall) metrics.

    Definition
    ----------
    This metric measures precision and recall on fingerprints distributions between
    reference and computed structures.

    - Precision (COV-P) measures the percentage of computed structures being close enough to
    the reference structures distribution. In other words, the number of computed structures
    that are similar to reference structures with respect to their fingerprints.

    - Recall (COV-R) measures the percentage of reference structures being close enough to
    the computed structures distribution. In other words, the number of reference structures
    that are similar to computed structures with respect to their fingerprints.

    Reference
    ---------
    Definition of Coverage metrics for crystals, using CrystalNN generated fingerprints:
    - Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
    Crystal Diffusion Variational Autoencoder for Periodic Material Generation.
    Bulletin of the American Physical Society, 67.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        ref_structs: list[Structure],
        transform: StructureDistribution,
        compute_precision: bool = True,
        compute_recall: bool = True,
        threshold: float = 0.4,
        **kwargs
    ) -> None:
        """
        Compute Coverage (Precision) and Coverage (Recall) metrics.
        
        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Coverage metrics on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        transform: StructureDistribution
            A callable taking a list of Structure objects and eventual keyword arguments and
            returning a numpy arrays of structures fingerprints.

        compute_precision: bool
            Whether to compute the Precision metric. Defaults to True.

        compute_recall: bool
            Whether to compute the Recall metric. Defaults to True.

        threshold: float
            Max distance to consider a data point as close enough to compared distribution.
            Defaults to 0.4.

        kwargs: Any
            Any additional keyword arguments to pass to the `transform` callable.
        """
        super().__init__(structures)
        self._ref_structs = ref_structs
        self._transform = ft.partial(transform, **kwargs)
        self._compute_precision = compute_precision
        self._compute_recall = compute_recall
        self._threshold = threshold

        self._compute()

    @staticmethod
    def _get_distance_closest(source: np.ndarray, target: np.ndarray) -> np.ndarray:
        """
        Computes the smallest distance between two structures distributions
        represented as numpy arrays.
        """
        src = torch.from_numpy(source)
        tgt = torch.from_numpy(target)

        idx_src, idx_tgt = knn(tgt, src, 1)
        closest_distance = (src[idx_src] - tgt[idx_tgt]).norm(dim=1)
        return closest_distance.numpy()

    def get_precision(self, source: np.ndarray, target: np.ndarray, threshold: float) -> float:
        """
        Get Precision metric value between two structures distributions
        represented as numpy arrays using initialized distance threshold.
        """
        distance = self._get_distance_closest(source=source, target=target)
        mask = distance < threshold
        return mask.astype(np.float32).mean().item()

    def get_recall(self, source: np.ndarray, target: np.ndarray, threshold: float) -> float:
        """
        Get Recall metric value between two structures distributions
        represented as numpy arrays using initialized distance threshold.
        """
        distance = self._get_distance_closest(source=target, target=source)
        mask = distance < threshold
        return mask.astype(np.float32).mean().item()


    def _get_metric_settings(self) -> dict[str, dict[str, tp.Any]]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "threshold": self._threshold,
                "transform": self._transform.func.__name__,
                "kwargs": self._transform.keywords,
                "compute_precision": self._compute_precision,
                "compute_recall": self._compute_recall
            }
        }

    def _compute(self) -> None:
        computed_fp = self._transform(self.structures)
        ref_fp = self._transform(self._ref_structs)

        if self._compute_precision:
            self._precision = self.get_precision(computed_fp, ref_fp, self._threshold)
        else:
            self._precision = None

        if self._compute_recall:
            self._recall = self.get_recall(computed_fp, ref_fp, self._threshold)
        else:
            self._recall = None

    @property
    def computed_data(self) -> list[MetricsData]:
        """
        List of all computed `MetricsData` objects containing structures from tested distribution
        and Coverage computation details (same values for all structures).
        """
        data = self._get_metric_settings()
        data[f"{type(self).__name__}_values"] = {
            "precision_percent": self.precision, "recall_percent": self.recall
        }
        return [
            MetricsData(
                struct,
                additional_data=data
            ) for struct in self.structures
        ]

    @property
    def precision(self) -> float | None:
        """Get computed Precision (COV-P) metric as a percentage."""
        return self._precision if self._precision is None else self._precision * 100

    @property
    def recall(self) -> float | None:
        """Get computed Recall (COV-R) metric as a percentage."""
        return self._recall if self._recall is None else self._recall * 100

    def write_result(self, filename: str, decimals: int = 6) -> None:
        if self._precision is None:
            precision_value = "Not computed"
        else:
            precision_value = f"{self.precision:.{decimals}f}%"

        if self._recall is None:
            recall_value = "Not computed"
        else:
            recall_value = f"{self.recall:.{decimals}f}%"

        text = "===== Coverage Results ====="
        text += f"Total computed structures:  {len(self)}"
        text += f"Total reference structures: {len(self._ref_structs)}"
        text += f"Precision (COV-P): {precision_value}"
        text += f"Recall (COV-R):    {recall_value}"

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))
