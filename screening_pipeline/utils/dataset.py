import json

from pymatgen.core import Element
from pymatgen.core.structure import Structure, SiteCollection, Composition


def load_phase_diagram_entries(filename: str) -> dict:
    with open(filename, "r") as fp:
        entries = json.load(fp)

    for entry in entries.values():
        entry["composition"] = Composition.from_dict(entry["composition"])

    return entries
