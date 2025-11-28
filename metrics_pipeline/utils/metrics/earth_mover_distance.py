"""Compute Earth Mover's Distance (or Wassertstein-1 Distance) metric."""

from pymatgen.core import Structure

from .metric_base import Metric


class EMD(Metric):
    """
    Compute Earth Mover's Distance (or Wassertstein-1 Distance) metric.
    
    Definition
    ----------
    TODO.

    For more details, see the original paper defining it:
    - Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
    Crystal Diffusion Variational Autoencoder for Periodic Material Generation.
    Bulletin of the American Physical Society, 67.
    """
    def __init__(self, structures: list[Structure], ref_structs: list[Structure], property: str) -> None:
        """
        Compute Earth Mover's Distance (or Wassertstein-1 Distance) metric.
        
        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Earth Mover's Distance metric on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        property: str
            Key inside structures properties to search for to compare them.
        """
        super().__init__(structures)
        self.ref_structs = ref_structs
        self.property = property

        self._compute()

    def _compute(self) -> None:
        raise NotImplementedError

    def write_result(self, filename: str, verbose: bool = False) -> None:
        return self._write_similarity_metric_result()