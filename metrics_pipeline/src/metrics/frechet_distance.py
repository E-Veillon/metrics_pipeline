"""Compute Fréchet Distance metric."""

import typing as tp
import functools as ft

import numpy as np

from pymatgen.core import Structure

from .metric_base import StructureDistribution, Metric, MetricsData
from src.utils import GenMatStructure


class FrechetDistance(Metric):
    """
    Compute Fréchet Distance metric.
    
    Definition
    ----------
    TODO.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        ref_structs: list[Structure],
        transform: StructureDistribution,
        **kwargs
    ) -> None:
        """
        Compute Fréchet Distance metric.
        
        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Fréchet Distance metric on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        transform: StructureDistribution
            Any callable taking a list of GenMatStructure objects and eventual keyword arguments
            and returning a numpy array representation of the structures distribution.
    
        kwargs: Any
            Any additional keyword arguments to pass to the `transform` callable.
        """
        super().__init__(structures)
        self.ref_structs = ref_structs
        self.transform = ft.partial(transform, **kwargs)
        self.kwargs = kwargs

        self._compute()

    @staticmethod
    def _sqrt_tr(x: np.ndarray) -> float:
        """Trace of the square root eigenvalues of a matrix."""
        return np.sum(np.sqrt(np.linalg.eigvals(x).real))

    def get_frechet_distance(self, x: np.ndarray, y: np.ndarray) -> float:
        """Frechet Distance between two matrices x and y."""
        mu_x = np.mean(x, axis=0)
        mu_y = np.mean(y, axis=0)

        sigma_x = np.cov(x.T)
        sigma_y = np.cov(y.T)

        fid: np.ndarray = (
            np.linalg.norm(mu_x - mu_y) ** 2
            + sigma_x.trace()
            + sigma_y.trace()
            - 2 * self._sqrt_tr(sigma_x @ sigma_y)
        )

        return fid.item()

    def _get_metric_settings(self) -> dict[str, dict[str, tp.Any]]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "transform": self.transform,
                "kwargs": self.kwargs
            }
        }

    def _compute(self) -> None:
        computed_array = self.transform(self.structures)
        ref_array = self.transform(self.ref_structs)
        self._distance = self.get_frechet_distance(computed_array, ref_array)

    @property
    def computed_data(self) -> list[MetricsData]:
        """
        List of all computed `MetricsData` objects containing structures from tested distribution
        and Fréchet Distance computation details (same values for all structures).
        """
        data = self._get_metric_settings()
        data[f"{type(self).__name__}_values"].update({"frechet_distance": self.computed_distance})
        return [
            MetricsData(
                struct,
                additional_data=data
            ) for struct in self.structures
        ]

    @property
    def computed_distance(self) -> float:
        """Get computed Fréchet Distance value."""
        return self._distance

    def write_result(self, filename: str, decimals: int = 6) -> None:
        text = "===== Fréchet Distance Results ====="
        text += f"Total computed structures:  {len(self)}"
        text += f"Total reference structures: {len(self.ref_structs)}"
        text += f"Fréchet Distance: {self.computed_distance:.{decimals}f}"

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))
