import json
from pymatgen.core.structure import Composition


def load_phase_diagram_entries(filename: str) -> dict:
    """Load entries data from a JSON file."""
    print(f"Loading {filename}...")
    with open(filename, "r") as fp:
        entries = json.load(fp)

    for entry in entries.values():
        entry["composition"] = Composition.from_dict(entry["composition"])
    print(f"Loading finished, {len(entries)} entries found.")
    return entries


