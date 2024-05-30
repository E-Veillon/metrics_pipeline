import json
from pymatgen.core.structure import Composition


def load_phase_diagram_entries(filename: str) -> dict:
    """Load entries data from a JSON file."""
    print("Loading reference dataset...")
    with open(filename, "r") as fp:
        entries = json.load(fp)

    for entry in entries.values():
        entry["composition"] = Composition.from_dict(entry["composition"])
    print("Loading finished.")
    return entries


