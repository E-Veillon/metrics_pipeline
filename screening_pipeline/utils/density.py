"""Functions to compute the density of Structure objects."""


from typing import List
from pymatgen.core import Structure
import numpy as np


def _volume_cm3(s: Structure) -> float:
    """Get the structure volume in cm^3."""
    return s.volume*1e-24

def _mass_g(s: Structure) -> float:
    """Get the structure mass in grams."""
    return sum(s.atomic_mass.to("g") for s in s.species)

def get_densities(structures: List[Structure]) -> np.ndarray:
    """Computes structure volumic mass in g / cm^3."""
    return np.array([_mass_g(s)/_volume_cm3(s) for s in structures])
