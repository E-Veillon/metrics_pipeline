"""Compute Fréchet ALIGNN Distance metric."""

from pymatgen.core import Structure

from .metric_base import Metric


class FAD(Metric):
    """
    Compute Fréchet ALIGNN Distance metric.
    
    Definition
    ----------
    TODO.

    Reference
    ---------
    Klipfel, A., Fregier, Y., Sayede, A., & Bouraoui, Z. (2024, March).
    Vector Field Oriented Diffusion Model for Crystal Material Generation.
    In Proceedings of the AAAI Conference on Artificial Intelligence (Vol. 38, No. 20, pp. 22193-22201).
    """
    def __init__(self, structures: list[Structure], ref_structs: list[Structure]) -> None:
        """
        Compute Fréchet ALIGNN Distance metric.
        
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