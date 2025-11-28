"""Compute Coverage (Precision) and Coverage (Recall) metrics."""

from pymatgen.core import Structure

from .metric_base import Metric


class Coverage(Metric):
    """
    Compute Coverage (Precision) and Coverage (Recall) metrics.
    
    Definition
    ----------
    This metric converts structures into fingerprints using the CrystalNN model, then compares
    distributions between the fingerprints of reference and computed structures.

    - Precision (COV-P) measures the proportion of computed structures being inside
    the reference structures distribution. In other words, the number of computed structures
    that are similar to reference structures with respect to their CrystalNN fingerprints.

    - Recall (COV-R) measures the proportion of reference structures being inside
    the computed structures distribution. In other words, the number of reference structures
    that are similar to computed structures with respect to their CrystalNN fingerprints.

    Reference
    ---------
    Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
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