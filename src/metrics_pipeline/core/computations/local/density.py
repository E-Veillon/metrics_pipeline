#!/usr/bin/python
"""Compute the density of Structure objects in g / cm³."""

import numpy as np

from pymatgen.core import Structure

from metrics_pipeline.core.utils.common_asserts import check_type


def _volume_cm3(structure: Structure) -> float:
    """Get the structure volume in cm³."""
    return structure.volume*1e-24


def _mass_g(structure: Structure) -> float:
    """Get the structure mass in grams."""
    return sum(s.atomic_mass.to("g") for s in structure.species)


def get_density(structure: Structure) -> float:
    """Get the structure volumic mass in g.cm⁻³."""
    return _mass_g(structure) / _volume_cm3(structure)


def get_densities(structures: list[Structure]) -> np.ndarray:
    """Computes structure volumic masses in g.cm⁻³."""
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))

    return np.array([get_density(s) for s in structures])
