from pymatgen.core import Structure
import numpy as np

from typing import List

def _volume_cm3(s:Structure)->float:
    return s.volume*1e-24

def _mass_g(s:Structure)->float:
    return sum([s.atomic_mass.to("g") for s in s.species])

def get_densities(
    structures: List[Structure]
) -> np.ndarray:
    return np.array([_mass_g(s)/_volume_cm3(s) for s in structures])

