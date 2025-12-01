"""Compute Coverage (Precision) and Coverage (Recall) metrics."""

import numpy as np
from tqdm.contrib.concurrent import process_map
from matminer.featurizers.site.fingerprint import CrystalNNFingerprint

import torch
from torch_cluster import knn

from pymatgen.core import Structure

from .metric_base import Metric


class Coverage(Metric):
    """
    Compute Coverage (Precision) and Coverage (Recall) metrics.
    
    Definition
    ----------
    This metric converts structures into fingerprints using the CrystalNN model, then compares
    distributions between the fingerprints of reference and computed structures.

    - Precision (COV-P) measures the proportion of computed structures being inside
    the reference structures distribution. In other words, the number of computed structures
    that are similar to reference structures with respect to their CrystalNN fingerprints.

    - Recall (COV-R) measures the proportion of reference structures being inside
    the computed structures distribution. In other words, the number of reference structures
    that are similar to computed structures with respect to their CrystalNN fingerprints.

    Reference
    ---------
    Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
    Crystal Diffusion Variational Autoencoder for Periodic Material Generation.
    Bulletin of the American Physical Society, 67.
    """
    _CrystalNNFP = CrystalNNFingerprint.from_preset("ops")

    def __init__(
        self,
        structures: list[Structure],
        ref_structs: list[Structure],
        threshold: float = 0.4,
        compute_precision: bool = True,
        compute_recall: bool = True,
        workers: int | None = None,
    ) -> None:
        """
        Compute Coverage (Precision) and Coverage (Recall) metrics.
        
        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Coverage metrics on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        threshold: float
            TODO. Defaults to 0.4.

        compute_precision: bool
            Whether to compute the Precision metric. Defaults to True.

        compute_recall: bool
            Whether to compute the Recall metric. Defaultsto True.

        workers: int, optional
            Number of parallel processes to spawn when computing fingerprints.
            If not given, `tqdm.contrib.concurrent.process_map()` default is used.
        """
        super().__init__(structures)
        self.ref_structs = ref_structs
        self.threshold = threshold
        self.compute_precision = compute_precision
        self.compute_recall = compute_recall
        self.workers = workers

        self._compute()

    def get_fingerprint(self, structure: Structure) -> np.ndarray | None:
        """Get the CrystalNN fingerprint of a structure."""
        try:
            atom_fingerprints = [
                self._CrystalNNFP.featurize(structure, i) for i, _ in enumerate(structure)
            ]
        except Exception:
            return None

        return np.mean(atom_fingerprints, axis=0)

    def to_crystalnn_fingerprints(self, structures: list[Structure]) -> list[np.ndarray | None]:
        """
        Convert given structures to their CrystalNN fingerprints.
        Can be parallelized over structures.

        Parameters
        ----------
        structures: list[Structure]
            Structures whoes fingerprint isneeded.

        Returns
        -------
        list[ndarray]
            List of CrystalNN fingerprints corresponding to given structures.
        """
        nb_structs = len(structures)
        chunksize = (min(nb_structs // 100, 10) if nb_structs >= 200 else 1)

        fingerprints = process_map(
            self.get_fingerprint,
            structures,
            max_workers=self.workers,
            chunksize=chunksize,
            desc="Convert structures to CrystalNN fingerprints"
        )
        return fingerprints

    @staticmethod
    def _get_distance_closest(source: np.ndarray, target: np.ndarray) -> np.ndarray:
        """Computes the smallest distance between two structures as coordinate arrays."""
        src = torch.from_numpy(source)
        tgt = torch.from_numpy(target)

        idx_src, idx_tgt = knn(tgt, src, 1)
        closest_distance = (src[idx_src] - tgt[idx_tgt]).norm(dim=1)
        return closest_distance.numpy()

    def get_precision(self, source: np.ndarray, target: np.ndarray, threshold: float) -> float:
        """Get Precision metric value."""
        distance = self._get_distance_closest(source=source, target=target)
        mask = distance < threshold
        return mask.astype(np.float32).mean().item()

    def get_recall(self, source: np.ndarray, target: np.ndarray, threshold: float) -> float:
        """Get Recall metric value."""
        distance = self._get_distance_closest(source=target, target=source)
        mask = distance < threshold
        return mask.astype(np.float32).mean().item()

    def _compute(self) -> None:
        computed_fp = self.to_crystalnn_fingerprints(self.structures)
        ref_fp = self.to_crystalnn_fingerprints(self.ref_structs)

        # TODO: explain why we eliminate a pair each time
        ref_fp, computed_fp = map(
            np.array,
            zip(
                *filter(
                    lambda x: x[0] is not None and x[1] is not None,
                    zip(ref_fp, computed_fp),
                )
            )
        )
        if self.compute_precision:
            self._precision = self.get_precision(computed_fp, ref_fp, self.threshold)
        else:
            self._precision = None

        if self.compute_recall:
            self._recall = self.get_recall(computed_fp, ref_fp, self.threshold)
        else:
            self._recall = None

    @property
    def precision(self) -> float | None:
        """Get computed Precision (COV-P) metric as a percentage."""
        return self._precision if self._precision is None else self._precision * 100

    @property
    def recall(self) -> float | None:
        """Get computed Recall (COV-R) metric as a percentage."""
        return self._recall if self._recall is None else self._recall * 100

    def write_result(self, filename: str, decimals: int = 6) -> None:
        if self._precision is None:
            precision_value = "Not computed"
        else:
            precision_value = f"{self.precision:.{decimals}f}%"

        if self._recall is None:
            recall_value = "Not computed"
        else:
            recall_value = f"{self.recall:.{decimals}f}%"

        text = "===== Coverage Results ====="
        text += f"Total computed structures:  {len(self)}"
        text += f"Total reference structures: {len(self.ref_structs)}"
        text += f"Precision (COV-P): {precision_value}"
        text += f"Recall (COV-R):    {recall_value}"

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))
