#!/usr/bin/python
"""Functions to compute vectors using a pretrained ALIGNN model."""


from typing import List, Literal, Optional
import numpy as np
import torch
from torch_geometric.data import Dataset
from materials_toolkit.data import StructureData, StructureLoader, collate
from materials_toolkit.models.alignn import get_pretrained_alignn
import tqdm

# PYTHON MATERIALS GENOMICS
from pymatgen.core import Structure, Element, Species

# LOCAL IMPORTS
from .common_asserts import check_type


########################################


def _species_to_tensor(elements: List[Element | Species]):
    """Convert Element objects to a Tensor containing their atomic numbers."""
    for idx, elt in enumerate(elements):
        if isinstance(elt, Species):
            elements[idx] = elt.element
    return torch.tensor([e.Z for e in elements], dtype=torch.long)


########################################


class StructuresDataset(Dataset):
    """
    A class to store several structures data
    and extract them into StructureData objects.
    """
    data_class = StructureData

    def __init__(self, structures: List[Structure]):
        super().__init__()

        data = [
            (
                torch.tensor(s.frac_coords, dtype=torch.float32),
                _species_to_tensor(s.species),
                torch.tensor(s.lattice.matrix.reshape(1, 3, 3), dtype=torch.float32),
            )
            for s in structures
        ]
        self.x = [x for x, _, _ in data]
        self.z = [z for _, z, _ in data]
        self.cell = [cell for _, _, cell in data]

    def len(self) -> int:
        """Number of structures in the instance."""
        return len(self.cell)

    def get(self, idx: int | torch.LongTensor) -> StructureData:
        """Extract structures data from indices."""
        if isinstance(idx, torch.LongTensor):
            return collate(
                [
                    StructureData(z=self.z[i], pos=self.x[i], cell=self.cell[i]) # type: ignore
                    for i in idx
                ]
            )
        return StructureData(z=self.z[idx], pos=self.x[idx], cell=self.cell[idx]) # type: ignore


########################################


@torch.no_grad
def vectors_from_alignn(
    structures: List[Structure],
    batch_size: int = 128,
    device: torch.device | str | None = None,
    model_name: str = "mp/e_form",
    output: Literal["latent","energy"] = "latent"
) -> np.ndarray:
    """
    Computes vector representation of structures with ALIGNN
    (latent or energy representation).
    """
    check_type(structures, "structures", (List,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))
    check_type(batch_size, "batch_size", (int,))
    if device is not None:
        check_type(device, "device", (torch.device,))
    check_type(model_name, "model_name", (str,))
    assert output in {"latent", "energy"}

    if device is None:
        if torch.cuda.is_available():
            device = "cuda"
        else:
            device = "cpu"

    alignn = get_pretrained_alignn(model_name).to(device) # type: ignore

    dataset = StructuresDataset(structures)
    loader = StructureLoader(dataset, batch_size=batch_size)
    latent = []

    def is_latent(output: str) -> bool:
        """True if 'output' == "latent"."""
        return output == "latent"

    for batch in tqdm.tqdm(loader):
        batch = batch.to(device)
        batch.build_graph(knn=12)
        batch.build_tripets()
        latent.append(alignn(batch, latent=is_latent(output)).detach())

    return torch.cat(latent).numpy()


########################################
