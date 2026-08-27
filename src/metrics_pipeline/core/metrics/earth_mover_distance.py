"""Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric."""

import typing as tp
import warnings
import functools as ft

import numpy as np
from scipy.stats import wasserstein_distance

from pymatgen.core import Structure

from .metric_base import StructureDistribution, Metric, MetricsData
from core.utils import GenMatStructure, StructureLike


class EMD(Metric):
    """
    Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric.

    Definition
    ----------
    Measure similarity between 1D distributions of one of the properties of
    reference and computed structures.

    Reference
    ---------
    Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
    Crystal Diffusion Variational Autoencoder for Periodic Material Generation.
    Bulletin of the American Physical Society, 67.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        ref_structs: list[Structure],
        transform: StructureDistribution | None = None,
        computed_property: str | None = None,
        **kwargs
    ) -> None:
        """
        Compute Earth Mover's Distance (or Wasserstein-1 Distance) metric.

        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Earth Mover's Distance metric on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        transform: StructureDistribution, optional
            Any callable taking a list of GenMatStructure objects and eventual keyword arguments
            and returning a numpy array representation of the structures 1D distribution
            of the computed property. If not given, `computed_property` must be given.

        computed_property: str, optional
            GenMatStructure property to compare. The property must be one of property attributes
            of GenMatStructure objects and must be defined for all computed and reference
            structures. If not given, `transform` must be given.

        kwargs: Any
            Additional keyword arguments to pass to `transform` (if used).

        Notes
        -----
        - If both `transform` and `computed_property` are given, `computed_property`
        is used in priority to save computations, and if it is not sufficient
        (i.e. not all structures have a valid value for the property), `transform`
        is used instead and stored property is ignored.
        """
        super().__init__(structures)

        def has_property(struct: StructureLike, property_name: str) -> bool:
            return (
                getattr(struct, property_name, None) is not None or
                struct.properties.get(property_name) is not None
            )

        if transform is None and computed_property is None:
            raise ValueError("Either 'transform' or 'computed_property' must be given.")

        if computed_property is not None:
            try:
                if not all(has_property(struct, computed_property) for struct in self.structures):
                    raise KeyError(
                        "Some computed structures do not have the property "
                        f"{computed_property!r} defined."
                    )
                if not all(has_property(struct, computed_property) for struct in ref_structs):
                    raise KeyError(
                        "Some reference structures do not have the property "
                        f"{computed_property!r} defined."
                    )
            except KeyError as exc:
                if transform is None:
                    raise exc
                else:
                    warn_msg = (
                        "Caught following exception while checking property "
                        f"'{computed_property}': {type(exc).__name__}: {exc} "
                        "Provided 'transform' callable will be used instead."
                    )
                    warnings.warn(warn_msg)
                    computed_property = None

        self._ref_structs = ref_structs
        self._transform = ft.partial(transform, **kwargs) if transform is not None else None
        self._computed_property = computed_property

        self._compute()

    def _get_metric_settings(self) -> dict[str, dict[str, tp.Any]]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "computed_property": self._computed_property,
                "transform": self._transform.func.__name__ if self._transform is not None else None,
                "kwargs": self._transform.keywords if self._transform is not None else None
            }
        }

    def _compute(self) -> None:
        if self._computed_property is None:
            if self._transform is None:
                raise RuntimeError("Type checker assertion.")
            computed_values = self._transform(self.structures)
            ref_values = self._transform(self._ref_structs)

        else:
            computed_values = np.array(
                [getattr(struct, self._computed_property) for struct in self.structures]
            )
            ref_values = np.array(
                [struct.properties[self._computed_property] for struct in self._ref_structs]
            )
        self._distance = wasserstein_distance(ref_values, computed_values)

    @property
    def computed_data(self) -> list[MetricsData]:
        """
        List of computed MetricsData objects containing structures from tested distribution
        and EMD computation details (same details for all structures).
        """
        data = self._get_metric_settings()
        data[f"{type(self).__name__}_values"] = {"distance": self.computed_distance}
        return [
            MetricsData(
                struct,
                additional_data=data
            ) for struct in self.structures
        ]

    @property
    def computed_distance(self) -> float:
        """Get computed EMD value."""
        return self._distance

    def write_result(self, filename: str, decimals: int = 6) -> None:
        text = "===== Earth Mover's Distance Results ====="
        text += f"Total computed structures:  {len(self)}"
        text += f"Total reference structures: {len(self._ref_structs)}"
        text += f"EMD ({self._computed_property}): {self.computed_distance:.{decimals}f}"

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))