#!/usr/bin/python
"""Functions to compute structures fingerprints using CrystalNN."""

from typing import List
from matminer.featurizers.site.fingerprint import CrystalNNFingerprint
import numpy as np
from tqdm.contrib.concurrent import process_map

# PYTHON MATERIALS GENOMICS
from pymatgen.core import Structure

# LOCAL IMPORTS
from .common_asserts import check_type, check_num_value


########################################


CrystalNNFP = CrystalNNFingerprint.from_preset("ops")

def _structure_to_fingerprint(struct: Structure):
    """Get the CrystalNN fingerprint of a structure."""
    try:
        atom_fingerprints = [
            CrystalNNFP.featurize(struct, i) for i, _ in enumerate(struct)
        ]
    except Exception:
        return None
    return np.mean(atom_fingerprints, axis=0)


########################################


def to_crystalnn_fingerprint(
    structures: List[Structure],
    workers: int = 1
) -> List[np.ndarray]:
    """
    Convert given structures to their CrystalNN fingerprints.
    Can be parallelized over structures.

    Parameters:
        structures ([Structure]):   Structures whoes fingerprint isneeded.

        workers (int):              Number of parallel processes to spawn.
    
    Returns: List[ndarray]:
        List of CrystalNN fingerprints corresponding to given structures.
    """
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))
    check_type(workers, "workers", (int,))
    check_num_value(workers, "workers", ">", 0)

    nb_structs = len(structures)
    chunksize = (min(nb_structs // 100, 10) if nb_structs >= 200 else 1)

    fingerprints = process_map(
        _structure_to_fingerprint,
        structures,
        max_workers=workers,
        chunksize=chunksize,
        desc="Convert structures to CrystalNN fingerprints"
    )
    return fingerprints


########################################
