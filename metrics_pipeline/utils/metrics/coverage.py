"""Compute Coverage (Precision) and Coverage (Recall) metrics."""

import functools as ft

import numpy as np
from tqdm.contrib.concurrent import process_map
from matminer.featurizers.site.fingerprint import CrystalNNFingerprint

import torch
from torch_cluster import knn

from pymatgen.core import Structure

from .metric_base import Metric, StructureFingerprint

# TODO: export to 'computations/models/crystalnn_fingerprints.py' when package is ready
def crystalnn_fingerprints(structures: list[Structure], workers: int | None = None):
    CrystalNNFP = CrystalNNFingerprint.from_preset("ops")
    def _get_fingerprint(structure: Structure) -> np.ndarray | None:
        try:
            atom_fingerprints = [
                CrystalNNFP.featurize(structure, i) for i, _ in enumerate(structure)
            ]
        except Exception:
            return None

        return np.mean(atom_fingerprints, axis=0)

    if workers is not None and workers == 0:
        fingerprints = [_get_fingerprint(struct) for struct in structures]
    else:
        nb_structs = len(structures)
        chunksize = (min(nb_structs // 100, 10) if nb_structs >= 200 else 1)

        fingerprints = process_map(
            _get_fingerprint,
            structures,
            max_workers=workers,
            chunksize=chunksize,
            desc="Convert structures to CrystalNN fingerprints"
        )

    return fingerprints

class Coverage(Metric):
    """
    Compute Coverage (Precision) and Coverage (Recall) metrics.

    Definition
    ----------
    This metric measures precision and recall on fingerprints distributions between
    reference and computed structures.

    - Precision (COV-P) measures the percentage of computed structures being close enough to
    the reference structures distribution. In other words, the number of computed structures
    that are similar to reference structures with respect to their fingerprints.

    - Recall (COV-R) measures the percentage of reference structures being close enough to
    the computed structures distribution. In other words, the number of reference structures
    that are similar to computed structures with respect to their fingerprints.

    Reference
    ---------
    Definition of Coverage metrics for crystals, using CrystalNN generated fingerprints:
    - Xie, T., Fu, X., Ganea, O., Barzilay, R., & Jaakkola, T. (2022).
    Crystal Diffusion Variational Autoencoder for Periodic Material Generation.
    Bulletin of the American Physical Society, 67.
    """
    def __init__(
        self,
        structures: list[Structure],
        ref_structs: list[Structure],
        transform: StructureFingerprint,
        compute_precision: bool = True,
        compute_recall: bool = True,
        threshold: float = 0.4,
        **kwargs
    ) -> None:
        """
        Compute Coverage (Precision) and Coverage (Recall) metrics.
        
        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Coverage metrics on.

        ref_structs: list[Structure]
            Known structures to use as reference distribution.

        transform: StructureFingerprint
            A callable taking a list of Structure objects and eventual keyword arguments and
            returning a list of numpy arrays representing structures fingerprints. Can return
            `None` for structures that could not be converted (e.g. weird unphysical structures).

        compute_precision: bool
            Whether to compute the Precision metric. Defaults to True.

        compute_recall: bool
            Whether to compute the Recall metric. Defaults to True.

        threshold: float
            Max distance to consider a data point as close enough to compared distribution.
            Defaults to 0.4.

        kwargs: Any
            Any additional keyword arguments to pass to the `transform` callable.
        """
        super().__init__(structures)
        self.ref_structs = ref_structs
        self.transform = ft.partial(transform, **kwargs)
        self.compute_precision = compute_precision
        self.compute_recall = compute_recall
        self.threshold = threshold

        self._compute()

    @staticmethod
    def _sanitize_fingerprints(
        computed_fp_list: list[np.ndarray | None], ref_fp_list: list[np.ndarray | None]
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Eliminate pairs of matching indices when at least one is `None` and convert
        valid data into single numpy arrays.
        """
        # TODO: explain why we should eliminate a pair each time
        # instead of simply eliminating None values in each list.
        computed_fp, ref_fp = map(
            np.array,
            zip(
                *filter(
                    lambda x: x[0] is not None and x[1] is not None,
                    zip(computed_fp_list, ref_fp_list),
                )
            )
        )
        return computed_fp, ref_fp

    @staticmethod
    def _get_distance_closest(source: np.ndarray, target: np.ndarray) -> np.ndarray:
        """
        Computes the smallest distance between two structures distributions
        represented as numpy arrays.
        """
        src = torch.from_numpy(source)
        tgt = torch.from_numpy(target)

        idx_src, idx_tgt = knn(tgt, src, 1)
        closest_distance = (src[idx_src] - tgt[idx_tgt]).norm(dim=1)
        return closest_distance.numpy()

    def get_precision(self, source: np.ndarray, target: np.ndarray, threshold: float) -> float:
        """
        Get Precision metric value between two structures distributions
        represented as numpy arrays using initialized distance threshold.
        """
        distance = self._get_distance_closest(source=source, target=target)
        mask = distance < threshold
        return mask.astype(np.float32).mean().item()

    def get_recall(self, source: np.ndarray, target: np.ndarray, threshold: float) -> float:
        """
        Get Recall metric value between two structures distributions
        represented as numpy arrays using initialized distance threshold.
        """
        distance = self._get_distance_closest(source=target, target=source)
        mask = distance < threshold
        return mask.astype(np.float32).mean().item()

    def _compute(self) -> None:
        computed_fp_list = self.transform(self.structures)
        ref_fp_list = self.transform(self.ref_structs)

        computed_fp, ref_fp = self._sanitize_fingerprints(computed_fp_list, ref_fp_list)

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
