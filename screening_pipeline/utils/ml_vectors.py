from pymatgen.core import Structure, Element
import numpy as np
import torch
from torch_geometric.data import Dataset
from materials_toolkit.data import StructureData, StructureLoader, collate
from materials_toolkit.models.alignn import get_pretrained_alignn
import tqdm

from typing import List,Literal


def _species_to_tensor(elements: List[Element]):
    return torch.tensor([e.Z for e in elements], dtype=torch.long)


class StructuresDataset(Dataset):
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
        return len(self.cell)

    def get(self, idx: int | torch.LongTensor) -> StructureData:
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
    structures: List[Structure],
    batch_size: int = 128,
    device: torch.device = None,
    model_name: str = "mp_e_form_alignn",
    output: Literal["latent","energy"]="latent"
) -> np.ndarray:
    assert output in ("latent","energy")

    if device is None:
        if torch.cuda.is_available():
            device = "cuda"
        else:
            device = "cpu"

    alignn = get_pretrained_alignn(model_name).to(device)

    dataset = StructuresDataset(structures)
    loader = StructureLoader(dataset, batch_size=batch_size)
    latent = []
    for batch in tqdm.tqdm(loader):
        batch = batch.to(device)
        batch.build_graph(knn=12)
        batch.build_tripets()
        latent.append(alignn(batch, latent=(output=="latent")).detach())

    return torch.cat(latent).numpy()
