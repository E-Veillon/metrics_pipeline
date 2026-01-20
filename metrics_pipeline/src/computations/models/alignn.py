"""Compute vectors using a pretrained ALIGNN model."""

import typing as tp

import tqdm

import numpy as np

import torch
from torch_geometric.data import Dataset

from materials_toolkit.data import StructureData, StructureLoader, collate
from materials_toolkit.models.alignn.pretrained import get_pretrained_alignn, models_name

from pymatgen.core import Structure, Element, Species

from src.utils import VisualIterator


def _species_to_tensor(elements: list[Element | Species]) -> torch.Tensor:
    """Convert Element objects to a Tensor containing their atomic numbers."""
    for idx, elt in enumerate(elements):
        if isinstance(elt, Species):
            elements[idx] = elt.element
    return torch.tensor([e.Z for e in elements], dtype=torch.long)


class StructuresDataset(Dataset):
    """
    A class to store several structures data
    and extract them into StructureData objects.
    """
    data_class = StructureData
    x: list[torch.FloatTensor]
    z: list[torch.LongTensor]
    cell: list[torch.FloatTensor]

    def __init__(self, structures: list[Structure]):
        super().__init__()

        data = [
            (
                torch.tensor(s.frac_coords, dtype=torch.float32),
                _species_to_tensor(s.species),
                torch.tensor(s.lattice.matrix.reshape(1, 3, 3), dtype=torch.float32),
            )
            for s in structures
        ]
        self.x = [tp.cast(torch.FloatTensor, x) for x, _, _ in data]
        self.z = [tp.cast(torch.LongTensor, z) for _, z, _ in data]
        self.cell = [tp.cast(torch.FloatTensor, cell) for _, _, cell in data]

    def len(self) -> int:
        """Number of structures in the instance."""
        return len(self.cell)

    def get(self, idx: int | torch.LongTensor) -> StructureData:
        """Extract structures data from indices."""
        if isinstance(idx, torch.LongTensor):
            return collate(
                [
                    StructureData(z=self.z[i], pos=self.x[i], cell=self.cell[i])
                    for i in idx
                ]
            )
        return StructureData(z=self.z[idx], pos=self.x[idx], cell=self.cell[idx])


@torch.no_grad
def vectors_from_alignn(
    structures: list[Structure],
    batch_size: int = 128,
    device: torch.device | str | None = None,
    model_name: models_name = "mp/e_form",
    output: tp.Literal["latent","energy"] = "latent",
    load_bar: tp.Literal["tqdm", "local", "quiet"] = "tqdm"
) -> np.ndarray:
    """
    Computes vector representation of structures with ALIGNN
    (latent or energy representation).
    """
    assert output in {"latent", "energy"}, ValueError(
        f"'output' only supports 'energy' and 'latent', got {output!r}."
    )
    if load_bar is not None:
        assert load_bar in {"tqdm", "local", "quiet"}, ValueError(
            f"'load_bar' only supports 'tqdm', 'local' and 'quiet', got {load_bar!r}."
        )

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    alignn = get_pretrained_alignn(model_name).to(device)

    dataset = StructuresDataset(structures)
    loader = StructureLoader(dataset, batch_size=batch_size)
    latent = []

    def is_latent(output: str) -> bool:
        """True if 'output' == "latent"."""
        return output == "latent"

    description = f"Comuting ALIGNN {output} values"

    if load_bar == "tqdm":
        loader = tqdm.tqdm(loader, desc=description)
    elif load_bar == "local":
        loader = VisualIterator(loader, desc=description, unit="computed", percent=True)

    for batch in loader:
        batch = batch.to(device)
        batch.build_graph(knn=12)
        batch.build_tripets()
        latent.append(alignn(batch, latent=is_latent(output)).detach())

    return torch.cat(latent).numpy()
