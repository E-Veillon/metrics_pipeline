"""Compute Validity metric."""

from pymatgen.core import Structure

from .metric_base import Metric

class Validity(Metric):
    """
    Compute Validity metric.
    
    Definition
    ----------
    A structure is valid if none of its atoms are closer than a threshold,
    generally 0.5 angstroms.
    """
    def __init__(self, structures: list[Structure], valid_tol: float = 0.5) -> None:
        """
        Compute Validity metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Validity on.

        valid_tol: float
            Validity threshold for atoms distance.
        """
        super().__init__(structures)
        self.valid_tol = valid_tol

        self._compute()

    def is_valid(self, structure: Structure) -> bool:
        """Whether a structure is valid with set tolerance."""
        return structure.is_valid(tol=self.valid_tol)

    def _compute(self) -> None:
        self._valid_structs = list(filter(self.is_valid, self.structures))
        self._invalid_structs = list(filter(lambda s: not self.is_valid(s), self.structures))

    @property
    def valid_structs(self) -> list[Structure]:
        """List of valid structures."""
        return self._valid_structs

    @property
    def invalid_structs(self) -> list[Structure]:
        """List of invalid structures."""
        return self._invalid_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        self._write_filter_metric_result(
            self._valid_structs, self._invalid_structs, filename, verbose
        )