"""Compute Root Mean Squared Displacement (RMSD) metric."""

import numpy as np

import torch
from torch import Tensor
from torch_scatter import scatter_mean

from pymatgen.core import Structure

from .metric_base import Metric


class RMSD(Metric):
    """
    Compute Root Mean Squared Displacement (RMSD) metric.
    
    Definition
    ----------
    Root mean squre distance between atomic positions before and after geometry optimization
    of a structure. The final result is the average RMSD over all structure pairs.
    """
    _offset_range = torch.arange(-5, 6, dtype=torch.float32)
    _offsets = torch.stack(
        torch.meshgrid(_offset_range, _offset_range, _offset_range, indexing="xy"), dim=3
    ).view(-1, 3)

    def __init__(self, structures: list[Structure], relaxed_structs: list[Structure]) -> None:
        """
        Compute Root Mean Squared Displacement (RMSD) metric.
        
        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Root Mean Squared Displacement metric on.

        relaxed_structs: list[Structure]
            The same structures in the same order as in `structures`,
            but after geometry optimization.
        """
        super().__init__(structures)
        assert len(self.structures) == len(relaxed_structs), (
            f"'structures' ({len(self.structures)}) and 'relaxed_structs' "
            f"({len(relaxed_structs)}) must have same length."
        )
        self.relaxed_structs = relaxed_structs

        self._compute()

    def _get_shortest_paths(
        self,
        x_src: torch.FloatTensor,
        cell_src: torch.FloatTensor,
        x_dst: torch.FloatTensor,
        num_atoms: torch.LongTensor,
    ) -> Tensor:

        idx = torch.arange(num_atoms.shape[0], dtype=torch.long, device=num_atoms.device)
        batch = idx.repeat_interleave(num_atoms)

        offset = self._offsets.clone().to(x_src.device)
        offset_euc = torch.einsum("ij,ljk->lik", offset, cell_src)

        paths = x_dst[:, None] + offset_euc[batch] - x_src[:, None]

        distance = paths.norm(dim=2)

        shortest_idx = distance.argmin(dim=1)
        idx = torch.arange(shortest_idx.shape[0], dtype=torch.long, device=idx.device)
        shortest_path = paths[idx, shortest_idx]

        return shortest_path

    @staticmethod
    def _to_euc(x: torch.FloatTensor, cell: torch.FloatTensor, batch: torch.LongTensor):
        return torch.einsum("ij,ijk->ik", x, cell[batch])

    @staticmethod
    def _to_inner(x: torch.FloatTensor, cell: torch.FloatTensor, batch: torch.LongTensor):
        return torch.einsum("ij,ijk->ik", x, cell[batch].inverse())

    @staticmethod
    def _center_around_zero(x: torch.FloatTensor) -> torch.FloatTensor:
        return (x + 0.5) % 1.0 + 0.5 # type: ignore

    @staticmethod
    def _polar(a: torch.FloatTensor) -> tuple[torch.FloatTensor, torch.FloatTensor]:
        w, s, vh = torch.linalg.svd(a)
        u = w @ vh
        return u, (vh.mT.conj() * s[:, None, :]) @ vh

    def rmsd(
        self,
        cell_src: torch.FloatTensor,
        x_src: torch.FloatTensor,
        cell_dst: torch.FloatTensor,
        x_dst: torch.FloatTensor,
        num_atoms: torch.LongTensor,
    ) -> torch.FloatTensor:
        """Compute RMSD between two structures."""
        batch_atoms = torch.arange(cell_src.shape[0], dtype=torch.long).repeat_interleave(
            num_atoms
        )

        _, cell_src = self._polar(cell_src)
        _, cell_dst = self._polar(cell_dst)
        x_src = self._center_around_zero(x_src)
        x_dst = self._center_around_zero(x_dst)

        x_src_euc = self._to_euc(x_src, cell_src, batch_atoms) # type: ignore
        x_dst_euc = self._to_euc(x_dst, cell_dst, batch_atoms) # type: ignore

        paths = _get_shortest_paths(x_src_euc, cell_src, x_dst_euc, num_atoms) # type: ignore

        avg_path = scatter_mean(paths, batch_atoms, dim=0, dim_size=num_atoms.shape[0])

        distance = (paths - avg_path[batch_atoms]).pow(2).sum(dim=1)

        return scatter_mean( # type: ignore
            distance, batch_atoms, dim=0, dim_size=num_atoms.shape[0]
        ).sqrt()
    # TODO: probablement le contenu de _compute()
    def rmsd_from_structures(
        self, struct1: list[Structure], struct2: list[Structure]
    ) -> np.ndarray:
        """
        Computes RMSD of each given pair of structure (pairing by index matching),
        and returns all results in a single array.
        """
        num_atoms = torch.tensor([len(s) for s in struct1], dtype=torch.long)

        assert (
            num_atoms == torch.tensor([len(s) for s in struct2], dtype=torch.long)
        ).all()

        x_src = torch.cat(
            [torch.tensor(s.frac_coords, dtype=torch.float32) for s in struct1], dim=0
        )
        cell_src = torch.tensor([s.lattice.matrix for s in struct1], dtype=torch.float32)

        x_dst = torch.cat(
            [torch.tensor(s.frac_coords, dtype=torch.float32) for s in struct2], dim=0
        )
        cell_dst = torch.tensor([s.lattice.matrix for s in struct2], dtype=torch.float32)

        return self.rmsd(cell_src, x_src, cell_dst, x_dst, num_atoms).numpy() # type: ignore

    def _compute(self) -> None:
        self._distance: float

    @property
    def get_rmsd(self) -> float:
        """Get computed RMSD value."""
        return self._distance

    def write_result(self, filename: str, decimals: int = 2, verbose: bool = False) -> None:
        text = "===== Root Mean Squared Displacement Results ====="
        text += f"Total structures:  {len(self)}"
        text += f"Average RMSD: {self._distance:.{decimals}f}"

        if verbose:
            text += "[Insert list of individual RMSD for all structures]" # TODO

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))