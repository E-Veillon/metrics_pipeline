#!/usr/bin/python
"""Functions to compute Root Mean Square Displacement between Structures."""


from typing import Tuple, List
import torch
#import torch.nn.functional as F
from torch_scatter import scatter_mean
import numpy as np

# PYHTON MATERIALS GENOMICS
from pymatgen.core import Structure


offset_range = torch.arange(-5, 6, dtype=torch.float32)
offsets = torch.stack(
    torch.meshgrid(offset_range, offset_range, offset_range, indexing="xy"), dim=3
).view(-1, 3)


def _get_shortest_paths(
    x_src: torch.FloatTensor,
    cell_src: torch.FloatTensor,
    x_dst: torch.FloatTensor,
    num_atoms: torch.LongTensor,
) -> torch.Tensor:

    idx = torch.arange(num_atoms.shape[0], dtype=torch.long, device=num_atoms.device)
    batch = idx.repeat_interleave(num_atoms)

    offset = offsets.clone().to(x_src.device)
    offset_euc = torch.einsum("ij,ljk->lik", offset, cell_src)

    paths = x_dst[:, None] + offset_euc[batch] - x_src[:, None]

    distance = paths.norm(dim=2)

    shortest_idx = distance.argmin(dim=1)
    idx = torch.arange(shortest_idx.shape[0], dtype=torch.long, device=idx.device)
    shortest_path = paths[idx, shortest_idx]

    return shortest_path


def _to_euc(x: torch.FloatTensor, cell: torch.FloatTensor, batch: torch.LongTensor):
    return torch.einsum("ij,ijk->ik", x, cell[batch])


def _to_inner(x: torch.FloatTensor, cell: torch.FloatTensor, batch: torch.LongTensor):
    return torch.einsum("ij,ijk->ik", x, cell[batch].inverse())


def _center_around_zero(x: torch.FloatTensor) -> torch.FloatTensor:
    return (x + 0.5) % 1.0 + 0.5 # type: ignore


def _polar(a: torch.FloatTensor) -> Tuple[torch.FloatTensor, torch.FloatTensor]:
    w, s, vh = torch.linalg.svd(a)
    u = w @ vh
    return u, (vh.mT.conj() * s[:, None, :]) @ vh


def rmsd(
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

    _, cell_src = _polar(cell_src)
    _, cell_dst = _polar(cell_dst)
    x_src = _center_around_zero(x_src)
    x_dst = _center_around_zero(x_dst)

    x_src_euc = _to_euc(x_src, cell_src, batch_atoms) # type: ignore
    x_dst_euc = _to_euc(x_dst, cell_dst, batch_atoms) # type: ignore

    paths = _get_shortest_paths(x_src_euc, cell_src, x_dst_euc, num_atoms) # type: ignore

    avg_path = scatter_mean(paths, batch_atoms, dim=0, dim_size=num_atoms.shape[0])

    distance = (paths - avg_path[batch_atoms]).pow(2).sum(dim=1)

    return scatter_mean( # type: ignore
        distance, batch_atoms, dim=0, dim_size=num_atoms.shape[0]
    ).sqrt()


def rmsd_from_structures(
    struct1: List[Structure], struct2: List[Structure]
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

    return rmsd(cell_src, x_src, cell_dst, x_dst, num_atoms).numpy() # type: ignore
