"""Compute Coverage (Precision) and Coverage (Recall) metrics."""

from pymatgen.core import Structure

from .metric_base import Metric


class Coverage(Metric):
    """
    Compute Coverage (Precision) and Coverage (Recall) metrics.
    
    Definition
    ----------
    TODO.

    For more details, see the original paper defining it:
    - Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
    Crystal Diffusion Variational Autoencoder for Periodic Material Generation.
    Bulletin of the American Physical Society, 67.
    """
    def __init__(self, structures: list[Structure], ref_structs: list[Structure]) -> None:
        """
        Compute Coverage (Precision) and Coverage (Recall) metrics.
        
        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Coverage metrics on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.
        """
        super().__init__(structures)
        self.ref_structs = ref_structs

        self._compute()

    def _compute(self) -> None:
        raise NotImplementedError

    def write_result(self, filename: str, verbose: bool = False) -> None:
        return self._write_similarity_metric_result()