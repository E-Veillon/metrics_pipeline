"""Subpackage for computations done directly inside GenMat."""

from .delta_sol import DSolCalcType, DSolFunc, DSolInput, DSolStructure, batch_get_dsol_band_gaps
from .density import get_densities
from .phase_diagrams import get_lacking_elts_entries

__all__ = [
    "DSolCalcType", "DSolFunc", "DSolInput", "DSolStructure", "batch_get_dsol_band_gaps",
    "get_densities",
    "get_lacking_elts_entries"
]