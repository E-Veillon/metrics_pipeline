"""Compute structures fingerprints using a pretrained CrystalNN model."""
import numpy as np
from tqdm.contrib.concurrent import process_map
from matminer.featurizers.site.fingerprint import CrystalNNFingerprint

from pymatgen.core import Structure


CrystalNNFP = CrystalNNFingerprint.from_preset("ops")


def _get_fingerprint(structure: Structure) -> np.ndarray | None:
    """Get CrystalNN fingerprint representation of a Structure."""
    try:
        atom_fingerprints = [
            CrystalNNFP.featurize(structure, i) for i, _ in enumerate(structure)
        ]
    except Exception:
        return None

    return np.mean(atom_fingerprints, axis=0)


def get_crystalnn_fingerprints(
    structures: list[Structure], workers: int | None = None
) -> np.ndarray:
    """
    Convert a list of structures to their CrystalNN fingerprints.
    
    Parameters
    ----------
    structures: list[Structure]
        Structure objects to convert.

    workers: int, optional
        Number of parallel processes to spawn. If not given,
        `tqdm.contrib.concurrent.process_map()` default is used.
        Pass 0 to disable `process_map()` and execute sequentially.

    Returns
    -------
    NDArray(dtype=np.float32)
        Numpy matrix of concatenated structures fingerprints. Rare occurrences
        of structures that cannot be converted are removed from the return.
    """
    if workers is not None and workers == 0:
        fingerprints = [_get_fingerprint(struct) for struct in structures]

    else:
        fingerprints = process_map(
            _get_fingerprint,
            structures,
            max_workers=workers,
            chunksize=min(10, len(structures) // 100 + 1),
            desc="Convert structures to CrystalNN fingerprints"
        )
    fingerprints = np.stack([f for f in fingerprints if f is not None], axis=0)
    return fingerprints
