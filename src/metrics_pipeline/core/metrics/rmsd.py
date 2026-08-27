"""Compute Root Mean Squared Displacement (RMSD) metric."""

from typing import Any

import numpy as np

import torch
from torch import Tensor
from torch_scatter import scatter_mean

from pymatgen.core import Structure

from .metric_base import Metric, MetricsData
from core.utils import GenMatStructure


class RMSD(Metric):
    """
    Compute Root Mean Squared Displacement (RMSD) metric.
    
    Definition
    ----------
    Root mean squre distance between atomic positions before and after geometry optimization
    of a structure. The final result is the average RMSD over all structure pairs.
    """
    _offset_range = torch.arange(-5, 6, dtype=torch.float32)
    _offsets_cpu = torch.stack(
        torch.meshgrid(_offset_range, _offset_range, _offset_range, indexing="xy"), dim=3
    ).view(-1, 3)
    _offsets_cache: dict[torch.device, torch.Tensor] = {}

    @classmethod
    def _get_offsets(cls, device: torch.device) -> torch.Tensor:
        """
        Get offset tensor on the specified device, caching to avoid redundant transfers.
        
        Parameters
        ----------
        device: torch.device
            Target device (CPU or GPU).
            
        Returns
        -------
        torch.Tensor
            Offset tensor on the specified device.
        """
        if device not in cls._offsets_cache:
            cls._offsets_cache[device] = cls._offsets_cpu.to(device)
        return cls._offsets_cache[device]

    def __init__(
        self,
        structures: list[GenMatStructure],
        relaxed_structs: list[Structure]
    ) -> None:
        """
        Compute Root Mean Squared Displacement (RMSD) metric.
        
        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Root Mean Squared Displacement metric on.

        relaxed_structs: list[Structure]
            The same structures in the same order as in `structures`,
            but after geometry optimization.
        """
        super().__init__(structures)
        if not len(self.structures) == len(relaxed_structs):
            raise ValueError(
                f"'structures' ({len(self.structures)}) and 'relaxed_structs' "
                f"({len(relaxed_structs)}) must have same length."
            )
        num_atoms = torch.tensor([len(struct) for struct in self.structures], dtype=torch.long)
        num_relaxed_atoms = torch.tensor(
                [len(struct) for struct in relaxed_structs], dtype=torch.long
            )
        if not (num_atoms == num_relaxed_atoms).all():
            ValueError(
                "Some structures do not have matching number of atoms before and after "
                "relaxation, please check that each matching index contain the same structure "
                "before and after relaxation."
            )
        self._num_atoms = num_atoms
        self._relaxed_structs = relaxed_structs

        self._compute()

    def _get_shortest_paths(
        self,
        x_src: Tensor,
        cell_src: Tensor,
        x_dst: Tensor,
    ) -> Tensor:
        """
        Compute shortest path between positions before and after optimization,
        taking account of periodic boundary conditions.
        
        Parameters
        ----------
        x_src: Tensor
            Positions before optimization.

        cell_src: Tensor
            Lattice vectors matrix of the structure before optimization.

        x_dst: Tensor
            Positions after optimization.

        Returns
        -------
        Tensor
            Shortest atoms trajectories between initial and optimized positions.
        """
        idx = torch.arange(self._num_atoms.shape[0], dtype=torch.long, device=self._num_atoms.device)
        batch = idx.repeat_interleave(self._num_atoms)

        offset = self._get_offsets(x_src.device)
        offset_euc = torch.einsum("ij,ljk->lik", offset, cell_src)

        paths = x_dst[:, None] + offset_euc[batch] - x_src[:, None]

        distance = paths.norm(dim=2)

        shortest_idx = distance.argmin(dim=1)
        idx = torch.arange(shortest_idx.shape[0], dtype=torch.long, device=idx.device)
        shortest_path = paths[idx, shortest_idx]

        return shortest_path

    @staticmethod
    def _to_euc(x: Tensor, cell: Tensor, batch: Tensor):
        return torch.einsum("ij,ijk->ik", x, cell[batch])

    @staticmethod
    def _center_around_zero(x: Tensor) -> Tensor:
        return (x + 0.5) % 1.0 + 0.5

    @staticmethod
    def _polar(a: Tensor) -> tuple[Tensor, Tensor]:
        w, s, vh = torch.linalg.svd(a)
        u = w @ vh
        return u, (vh.mT.conj() * s[:, None, :]) @ vh

    def compute_rmsd(
        self,
        x_src: Tensor,
        cell_src: Tensor,
        x_dst: Tensor,
        cell_dst: Tensor,
    ) -> Tensor:
        """
        Compute RMSD between two structures represented as Tensors.
        
        Parameters
        ----------
        x_src: Tensor
            Fractional positions before optimization.

        cell_src: Tensor
            Lattice vectors matrix of the non-optimized structure.

        x_dst: Tensor
            Fractional positions after optimization.

        cell_dst: Tensor
            Lattice vectors matrix of the optimized structure.
        """
        batch_atoms = torch.arange(cell_src.shape[0], dtype=torch.long).repeat_interleave(
            self._num_atoms
        )
        _, cell_src = self._polar(cell_src)
        _, cell_dst = self._polar(cell_dst)
        x_src = self._center_around_zero(x_src)
        x_dst = self._center_around_zero(x_dst)

        x_src_euc = self._to_euc(x_src, cell_src, batch_atoms)
        x_dst_euc = self._to_euc(x_dst, cell_dst, batch_atoms)

        paths = self._get_shortest_paths(x_src_euc, cell_src, x_dst_euc)

        avg_path = scatter_mean(paths, batch_atoms, dim=0, dim_size=self._num_atoms.shape[0])

        distance = (paths - avg_path[batch_atoms]).pow(2).sum(dim=1)

        return scatter_mean(
            distance, batch_atoms, dim=0, dim_size=self._num_atoms.shape[0]
        ).sqrt()

    def _get_metric_settings(self) -> dict[str, Any]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {}
        }

    def _compute(self) -> None:
        x_src = torch.cat(
            [torch.tensor(struct.frac_coords, dtype=torch.float32) for struct in self.structures],
            dim=0
        )
        cell_src = torch.tensor(
            [struct.lattice.matrix for struct in self.structures], dtype=torch.float32
        )
        x_dst = torch.cat(
            [torch.tensor(struct.frac_coords, dtype=torch.float32) for struct in self._relaxed_structs],
            dim=0
        )
        cell_dst = torch.tensor(
            [struct.lattice.matrix for struct in self._relaxed_structs], dtype=torch.float32
        )
        self._distances = self.compute_rmsd(x_src, cell_src, x_dst, cell_dst).numpy(force=True)

    @property
    def computed_data(self) -> list[MetricsData]:
        """
        List of computed `MetricsData` objects each containing an unrelaxed structure,
        its associated RMSD value, and its relaxed version as additional data.
        """
        return [
            MetricsData(struct, rmsd=rmsd_val, additional_data={"relaxed_structure": relaxed})
            for struct, relaxed, rmsd_val in zip(
                self.structures, self._relaxed_structs, self.all_rmsd
            )
        ]

    @property
    def average_rmsd(self) -> float:
        """Get computed average RMSD value over all structures."""
        return float(np.mean(self._distances).item())

    @property
    def all_rmsd(self) -> np.ndarray:
        """Get 1D array of computed RMSD values for all structures."""
        return self._distances

    def get_rmsd_from_index(self, idx: int) -> float:
        """Get computed RMSD value of a specific structure from its index."""
        return float(self.all_rmsd[idx])

    def write_result(self, filename: str, decimals: int = 6, verbose: bool = False) -> None:
        text = "===== Root Mean Squared Displacement Results ====="
        text += f"Total structures:  {len(self)}"
        text += f"Average RMSD: {self.average_rmsd:.{decimals}f}"

        if verbose:
            max_digits = len(str(len(self.all_rmsd)))
            rmsd_text = "\n".join(
                [
                    f"{idx:0>{max_digits}} - {rmsd:.{decimals}f}"
                    for idx, rmsd in enumerate(self.all_rmsd)
                ]
            )
            text += f"List of RMSD for each structure:\n{rmsd_text}\n"

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))
