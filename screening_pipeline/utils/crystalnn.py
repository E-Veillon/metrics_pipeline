from typing import List
from pymatgen.core import Structure
from matminer.featurizers.site.fingerprint import CrystalNNFingerprint
import numpy as np
from tqdm.contrib.concurrent import process_map

CrystalNNFP = CrystalNNFingerprint.from_preset("ops")


def _structure_to_fingerprint(struct: Structure):
    """Get the CrystalNN fingerprint of a structure."""
    try:
        atom_fingerprints = [
            CrystalNNFP.featurize(struct, i) for i, _ in enumerate(struct)
        ]
    except:
        return None
    return np.mean(atom_fingerprints, axis=0)


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
    fingerprints = process_map(
        _structure_to_fingerprint,
        structures,
        max_workers=workers,
        desc="Convert structures to CrystalNN fingerprints"
    )
    return fingerprints
