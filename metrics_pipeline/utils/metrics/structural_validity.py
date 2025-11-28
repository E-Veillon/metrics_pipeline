"""Compute Structural Validity metric."""

from collections import OrderedDict

from pymatgen.core import Structure

from .metric_base import Metric

class StructValidity(Metric):
    """
    Compute Structural Validity metric.
    
    Definition
    ----------
    A structure is valid if none of its atoms are closer than a threshold, generally 0.5 angstroms.

    Reference
    ---------
    Court, C. J., Yildirim, B., Jain, A., & Cole, J. M. (2020).
    3-D inorganic crystal structure generation and property prediction via representation learning.
    Journal of Chemical Information and Modeling, 60(10), 4518-4535.
    """
    def __init__(self, structures: list[Structure], valid_tol: float = 0.5) -> None:
        """
        Compute Structural Validity metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Structural Validity on.

        valid_tol: float
           Structural Validity threshold for atoms distance.
           Defaults to 0.5 angstroms.
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
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Valid", self._valid_structs),
                ("Invalid", self._invalid_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)