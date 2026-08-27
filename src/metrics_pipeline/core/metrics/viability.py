"""Compute Viability metric."""

import itertools as itt
from collections import OrderedDict
from typing import Any

import numpy as np
from scipy.spatial.distance import pdist

from .metric_base import Metric, MetricsData
from .backend import (
    slater_radii_table_pm_1,
    clementi_et_al_radii_table_pm
)
from metrics_pipeline.core.utils.periodic_table import ALL_ELT_SYMBOL_TO_Z
from metrics_pipeline.core.utils.genmat_data import GenMatStructure, StructureLike
from metrics_pipeline.core.genmat_io.json import JsonLoader


class BadTableError(ValueError):
    """Loaded data is not a valid table."""
    ...


class Viability(Metric):
    """
    Compute Viability metric.

    Definition
    ----------
    A structure is viable if its lattice lengths are at least as long as the diameter of the
    biggest atom in the cell, its lattice angles are between 10° and 170°, and none of its
    atoms are closer than a threshold, based on a table of atomic radii for all 118 elements.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        table: str = "clementi",
        abs_min_dist: float = 0.5,
        workers: int | None = None
    ) -> None:
        """
        Compute Viability metric.

        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Viability on.

        table: str
            Name of the atomic radii table to use as reference.
            Supports:
            - 'slater' for experimental radii measured by J. C. Slater (1964).
            - 'clementi' (default) for simulated radii calculated by Clementi et al. (1963 and 1967).
            - Any file path to provide your own JSON file containing a dict of {'symbol': radius}
            with all 118 elements radii as integers in picometers.

        abs_min_dist: float
            Absolute minimal distance in angstroms to consider two atoms as too close,
            whatever their radii. Can be useful for unknown radii having a default value,
            to ensure at least a physical minimum distance comparison. Defaults to 0.5.

        workers: int, optional
            Number of parallel processes to spawn for high throughput structure matching.
            If not given, will use default of `tqdm.contrib.concurrent.process_map()`.
            pass 0 to disable `process_map()` and do it sequentially.
        """
        super().__init__(structures, workers)
        self._radii_table = self._get_radii_table(table)
        self._table_name = table
        self._abs_min_dist = abs_min_dist

        self._compute()

    @staticmethod
    def _get_radii_table(table_name: str) -> dict[str, float]:
        """Load a reference atomic radii table."""
        match table_name:
            case "slater":
                table = slater_radii_table_pm_1
            case "clementi":
                table = clementi_et_al_radii_table_pm
            case str():
                table = JsonLoader(table_name).load_as_dict()
            case _:
                raise TypeError(
                    f"'table_name' expected a type 'str', got {type(table_name).__name__!r}."
                )

        all_elts = set(ALL_ELT_SYMBOL_TO_Z.keys())

        if not all(isinstance(key, str) for key in table.keys()):
            key_types = ", ".join(sorted(type(key).__name__ for key in table.keys()))
            raise TypeError(
                f"Some keys in loaded table are not of type 'str', got {key_types}."
            )
        if not all(isinstance(val, int) for val in table.values()):
            val_types = ", ".join(sorted(type(val).__name__ for val in table.values()))
            raise TypeError(
                f"Some values in loaded table are not of type 'int', got {val_types}."
            )
        if not all(table.get(elt) is not None for elt in all_elts):
            lacking_elts = ", ".join(sorted(elt for elt in all_elts if table.get(elt) is None))
            raise BadTableError(
                f"Some elements were not found in loaded table: {lacking_elts}."
            )
        return {k: v / 100 for k, v in table.items()}

    def _get_min_dist(self, radius1: float, radius2: float) -> float:
        """Compute minimal distance between two atoms with respect to their radii."""
        min_radii_dist = sum((radius1, radius2)) * 0.9 # 10% uncertainty tolerance
        return max(self._abs_min_dist, round(min_radii_dist, 3)) # round result to 0.1 pm scale

    def is_viable(self, structure: StructureLike) -> bool:
        """
        Compute Viability of a structure. Compatible with standard `Structure` objects,
        but returns a simple bool instead of a computed `MetricsData` object.
        """
        if any(not (10 <= angle <= 170) for angle in structure.lattice.angles):
            return False

        radii = [self._radii_table[elt.symbol] for elt in structure.species]

        if any(length < 2 * max(radii) for length in structure.lattice.lengths):
            return False

        min_dists = np.array(list(itt.starmap(self._get_min_dist, itt.combinations(radii, r=2))))
        # Use pdist to compute only upper-triangle distances, avoiding redundant computation
        true_dists = pdist(structure.cart_coords)
        return np.all(np.subtract(true_dists, min_dists) >= 0.0).item()

    def get_viability(self, structure: GenMatStructure) -> MetricsData:
        """
        Compute Viability of a structure. Return a computed `MetricsData` containing
        the structure and metric results.
        """
        return MetricsData(
            structure,
            is_viable=self.is_viable(structure),
            additional_data=self._get_metric_settings()
        )

    def _get_metric_settings(self) -> dict[str, Any]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "abs_min_dist": self._abs_min_dist,
                "table_name": self._table_name
            }
        }

    def _compute(self) -> None:
        if self.workers == 0:
            self._computed_data = list(self._sequential_compute(
                    self.get_viability, self.structures
                ))
        else:
            self._computed_data = self._parallel_compute(
                self.get_viability, self.structures
            )

    @property
    def computed_data(self) -> list[MetricsData]:
        """List of all computed MetricsData objects."""
        return self._computed_data

    @property
    def viable_data(self) -> list[MetricsData]:
        """List of computed MetricsData containing a viable structure."""
        return [data for data in self.computed_data if data.is_viable]

    @property
    def non_viable_data(self) -> list[MetricsData]:
        """List of computed MetricsData containing a non-viable structure."""
        return [data for data in self.computed_data if data.is_viable is False]

    @property
    def viable_structs(self) -> list[GenMatStructure]:
        """List of viable structures."""
        return [data.typed_structure for data in self.computed_data if data.is_viable]

    @property
    def non_viable_structs(self) -> list[GenMatStructure]:
        """List of non-viable structures."""
        return [data.typed_structure for data in self.computed_data if data.is_viable is False]

    @property
    def viable_names(self) -> list[str]:
        """List of names of viable structures."""
        return [data.typed_structure.name for data in self.computed_data if data.is_viable]

    @property
    def non_viable_names(self) -> list[str]:
        """List of names of non-viable structures."""
        return [data.typed_structure.name for data in self.computed_data if data.is_viable is False]

    @property
    def viable_names_set(self) -> set[str]:
        """Set of names of viable structures."""
        return {data.typed_structure.name for data in self.computed_data if data.is_viable}

    @property
    def non_viable_names_set(self) -> set[str]:
        """Set of names of non-viable structures."""
        return {data.typed_structure.name for data in self.computed_data if data.is_viable is False}

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[GenMatStructure]] = OrderedDict(
            [
                ("Viable", self.viable_structs),
                ("Non viable", self.non_viable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
