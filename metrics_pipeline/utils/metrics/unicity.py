"""Compute Unicity metric."""

import itertools as itt

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from .metric_base import Metric


class Unicity(Metric):
    """
    Compute Unicity metric.
    
    Definition
    ----------
    Structures are uniques if they are not an equivalent representation of a previous structure
    in given list.
    """
    def __init__(
        self,
        structures: list[Structure],
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5.0
    ) -> None:
        """
        Compute Unicity metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Unicity on.

        ltol: float
            Fractional length tolerance for structures matching. Default is 0.2.

        stol: float
            Site tolerance for structures matching. Defined as the fraction of the
            average free length per atom, i.e. ( V / Nsites ) ** (1/3). Default is 0.3.

        angle_tol: float
            Angle tolerance for structures matching in degrees. Default is 5.0 degrees.
        """
        super().__init__(structures)
        self.matcher = StructureMatcher(
            ltol=ltol, stol=stol, angle_tol=angle_tol, scale=False, attempt_supercell=True
        )

        self._compute()

    def is_unique(self, structure: Structure) -> bool:
        """Whether given structure is equivalent to any stored unique structure."""
        return any(self.matcher.fit(structure, struct) for struct in self.unique_structs)

    def _compute(self) -> None:
        groups: list[list[Structure]] = self.matcher.group_structures(self.structures)
        self._unique_structs = [group[0] for group in groups]
        self._duplicate_structs = list(itt.chain.from_iterable([group[1:] for group in groups]))

    @property
    def unique_structs(self) -> list[Structure]:
        return self._unique_structs

    @property
    def duplicate_structs(self) -> list[Structure]:
        return self._duplicate_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        return self._write_filter_metric_result(
            self._unique_structs, self._duplicate_structs, filename, verbose
        )
        