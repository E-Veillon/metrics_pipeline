"""Compute Symmetry metric."""

import functools as ft
from collections import OrderedDict

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from .metric_base import Metric


class Symmetry(Metric):
    """
    Compute Symmetry metric.

    Definition
    ----------
    Structures are considered symmetric if their spacegroup symmetry is not one of the
    triclinic crystal systems, i.e. P1 (no space symmetry element) or P-1 (inversion center only).
    """
    def __init__(
        self,
        structures: list[Structure],
        symprec: float = 0.01,
        angleprec: float = 5.0
    ) -> None:
        """
        Compute Symmetry metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Symmetry on.

        symprec: float
            Tolerance for symmetry finding. See `pymatgen.symmetry.analyzer.SpacegroupAnalyzer`
            for more details. Defaults to 0.01.

        angleprec: float
            Angle tolerance for symmetry finding. See `pymatgen.symmetry.analyzer.SpacegroupAnalyzer`
            for more details. Defaults to 5.0 degrees.
        """
        super().__init__(structures)
        self.symprec = symprec
        self.angleprec = angleprec
        self.analyzer = ft.partial(SpacegroupAnalyzer, symprec=symprec, angle_tolerance=angleprec)

        self._compute()

    def is_symmetric(self, structure: Structure) -> bool:
        return self.analyzer(structure).get_crystal_system() != "triclinic"

    def _compute(self) -> None:
        self._symmetric_structs = [
            struct for struct in self.structures if self.is_symmetric(struct)
        ]
        self._triclinic_structs = [
            struct for struct in self.structures if not self.is_symmetric(struct)
        ]

    @property
    def symmetric_structs(self) -> list[Structure]:
        """List of structures of higher symmetry than a triclinic system."""
        return self._symmetric_structs

    @property
    def triclinic_structs(self) -> list[Structure]:
        """List of structures of triclinic crystal system."""
        return self._triclinic_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Symmetric", self._symmetric_structs),
                ("Triclinic", self._triclinic_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
