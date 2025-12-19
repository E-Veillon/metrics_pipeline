#!/usr/bin/python
"""Compute the density of Structure objects in g / cm³."""

import numpy as np

from pymatgen.core import Structure

from src.utils import check_type


def _volume_cm3(s: Structure) -> float:
    """Get the structure volume in cm³."""
    return s.volume*1e-24


def _mass_g(s: Structure) -> float:
    """Get the structure mass in grams."""
    return sum(s.atomic_mass.to("g") for s in s.species)


def get_densities(structures: list[Structure]) -> np.ndarray:
    """Computes structure volumic mass in g / cm³."""
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))

    return np.array([_mass_g(s)/_volume_cm3(s) for s in structures])
