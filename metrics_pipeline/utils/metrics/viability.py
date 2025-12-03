"""Compute Viability metric."""

import os
import json
import itertools as itt
from collections import OrderedDict

import numpy as np
from pymatgen.core import Structure, Element

from .metric_base import Metric
from .radii_table import (
    slater_radii_table_pm_1,
    clementi_et_al_radii_table_pm
)

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
        structures: list[Structure],
        table: str = "clementi",
        min_dist: float = 0.5
    ) -> None:
        """
        Compute Viability metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Viability on.

        table: str
            Name of the atomic radii table to use as reference.
            Supports:
            - 'slater' for experimental radii measured by J. C. Slater (1964).
            - 'clementi' (default) for simulated radii calculated by Clementi et al. (1963 and 1967).
            - Any file path to provide your own JSON file containing a dict of {'symbol': radius}
            with all 118 elements radii as integers in picometers.

        min_dist: float
            Absolute minimal distance in angstroms to consider two atoms as too close,
            whatever their radii. Can be useful for unknown radii having a default value,
            to ensure at least a physical minimum distance comparison. Defaults to 0.5.
        """
        super().__init__(structures)
        self.radii_table = self._get_radii_table(table)
        self.min_dist = min_dist

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
                assert os.path.isfile(table_name), FileNotFoundError(
                    f"{table_name}: No such file found."
                )
                with open(table_name, "rt", encoding="utf-8") as fp:
                    table = json.load(fp)
            case _:
                raise TypeError(
                    f"'table_name' expected a type 'str', got {type(table_name).__name__!r}."
                )

        all_elts = set(Element.__members__)
        if "D" in all_elts:
            all_elts.remove("D") # Remove Deuterium (H isotope)
        if "T" in all_elts:
            all_elts.remove("T") # Remove Tritium (H isotope)

        assert isinstance(table, dict)
        assert all(isinstance(key, str) for key in table.keys())
        assert all(isinstance(val, int) for val in table.values())
        assert all(table.get(elt) is not None for elt in all_elts)

        return {k: v / 100 for k, v in table.items()}

    def _get_min_dist(self, radius1: float, radius2: float) -> float:
        """Compute minimal distance between two atoms with respect to their radii."""
        min_radii_dist = sum((radius1, radius2)) * 0.9 # 10% uncertainty tolerance
        return max(self.min_dist, round(min_radii_dist, 3)) # round result to 0.1 pm scale

    def is_viable(self, structure: Structure) -> bool:
        if any(not (10 <= angle <= 170) for angle in structure.lattice.angles):
            return False

        radii = [self.radii_table[elt.symbol] for elt in structure.species]

        if any(length < 2 * max(radii) for length in structure.lattice.lengths):
            return False

        min_dists = np.array(list(map(self._get_min_dist, *itt.combinations(radii, r=2))))
        true_dists = structure.distance_matrix[np.triu_indices(len(structure), 1)]
        return np.all(np.subtract(true_dists, min_dists) >= 0.0).item()

    def _compute(self) -> None:
        self._viable_structs = [struct for struct in self.structures if self.is_viable(struct)]
        self._non_viable_structs = [struct for struct in self.structures if not self.is_viable(struct)]

    @property
    def viable_structs(self) -> list[Structure]:
        """List of viable structures."""
        return self._viable_structs

    @property
    def non_viable_structs(self) -> list[Structure]:
        """List of non-viable structures."""
        return self._non_viable_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Viable", self._viable_structs),
                ("Non viable", self._non_viable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
