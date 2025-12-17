"""Subpackage for computations done directly inside GenMat."""

from .density import get_densities
from .phase_diagrams import get_lacking_elts_entries

__all__ = [
    "get_densities", "get_lacking_elts_entries"
]