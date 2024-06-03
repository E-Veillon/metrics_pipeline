from typing import List
from pymatgen.core import Structure
from matminer.featurizers.site.fingerprint import CrystalNNFingerprint
import numpy as np

CrystalNNFP = CrystalNNFingerprint.from_preset("ops")


def _structure_to_fingerprint(struct: Structure):
    try:
        atom_fingerprints = [
            CrystalNNFP.featurize(struct, i) for i, _ in enumerate(struct)
        ]
    except:
        return None
    return np.mean(atom_fingerprints, axis=0)


def to_crystalnn_fingerprint(
    structures: List[Structure],
) -> List[np.ndarray]:
    return [_structure_to_fingerprint(structure) for structure in structures]
