"""Downloads and manages remote databases API and data."""

import os
import json
from typing import List, Dict, Any
from mp_api.client import MPRester
from emmet.core.thermo import ThermoType
from pymatgen.core import SETTINGS
from pymatgen.core.structure import Composition


########################################


class APIKeyNotFoundError(Exception):
    """
    An Error ocurring when an API key is not defined while trying to access an API.
    """


########################################


def _mp_api_key_check(api_key: str|None = None) -> str:
    """
    Check whether the MP API key is defined either in .pmgrc.yaml or manually.
    If it is found, the key string is returned. Otherwise, an exception is raised.
    """
    api_key = api_key or SETTINGS.get("PMG_MAPI_KEY")
    if not api_key:
        raise APIKeyNotFoundError(
            "An API key must be given to access MP API. "
            "Configure it in .pmgrc.yaml in 'PMG_MAPI_KEY' keyword or pass it manually."
        )
    return api_key


########################################


def mp_api_download(path: str, api_key: str|None = None) -> None:
    """
    Download GGA/GGA+U phase diagram relevant entries from MP API to a local JSON file.
    Downloaded entries are dicts of the form {"entry_id": str, "composition": dict, "final_energy": float}
    
    Parameters:
        path (str):     Where to build the output JSON file.

        api_key (str):  User's API KEY to access the Materials Project API data.
                        If not given directly, it is searched in .pmgrc.yaml with
                        the "PMG_MAPI_KEY" variable.
    """
    api_key = _mp_api_key_check(api_key)

    with MPRester(
        api_key=api_key,
        use_document_model=False
        ) as mpr:
        data = mpr.materials.thermo.search(thermo_types=ThermoType.GGA_GGA_U,
            all_fields=False, fields=["material_id","composition","energy_per_atom"]
        )
        print(f"Number of Materials Project entries downloaded: {len(data)}")

        final_data = {}
        for entry in data:
            nb_atoms = Composition(entry["composition"]).num_atoms
            final_data.update(
                {
                    entry["material_id"]: {
                        "entry_id": entry["material_id"],
                        "composition": entry["composition"],
                        "final_energy": entry["energy_per_atom"] * nb_atoms
                    }
                }
            )

    with open(path,"wt") as fp:
        json.dump(data, fp)


########################################


def mp_api_print_avail_fields(endpoint: str, api_key: str|None = None) -> None:
    """Print the list of available fields corresponding to given MPRester endpoint."""
    if not isinstance(endpoint, str):
        raise TypeError(
            f"'endpoint' arg expected a type 'str', got '{type(endpoint)}' instead."
        )

    api_key = _mp_api_key_check(api_key)

    with MPRester("ZIQWsd2I9w337cKTrEKbJA0XKcveFIk3") as mpr:
        if endpoint == "materials":
            print(mpr.materials.thermo.available_fields)
        elif endpoint == "thermo":
            print(mpr.materials.available_fields)
        else:
            raise NotImplementedError(
                "Given 'endpoint' arg does not match any implemented doc."
            )


########################################


def process_oqmd_json_file(path: str) -> None:
    """
    Process OQMD json data file to get needed data format in a new JSON file.
    Expected OQMD JSON key formatting is the one got with OQMD RESTful API.
    The processed file have a '_proc' suffix added to its name.

    Parameters:
        path (str): Path to the OQMD JSON file to process.
    """
    assert os.path.isfile(path), f"{path}: No such file found."
    assert path.endswith(".json"), (
        f"Given filename '{os.path.basename(path)}' is not a valid JSON file."
    )

    with open(path, "rt") as fp:
        data = json.load(fp)
    
    processed_data = {}
    
    for name, struct in data.items():
        struct: Dict[str, Any]

        if (
            struct.get("entry_id") is None
            or struct.get("composition") is None
            or struct.get("natoms") is None
            or struct.get("delta_e") is None
        ):
            if name.startswith("proc-oqmd-"):
                raise ValueError(
                    f"structure '{name}' was already processed (begins with 'proc-oqmd-'). "
                    "Please verify that you are not processing the same data more than "
                    "once, as it may cause wrong total energy computations."
                )
            else:
                raise ValueError(
                    f"At least one of needed data keys is missing in structure '{name}'. "
                    "Please verify that the following keys are present: "
                    "'entry_id', 'composition', 'natoms', 'delta_e'."
                )

        entry_id: str = struct.get("entry_id")

        if isinstance(struct.get("composition"), Dict):
            composition: Dict[str, int] = struct.get("composition")

        elif isinstance(struct.get("composition"), str):
            old_comp = Composition("".join(struct.get("composition").split()), strict=True)

            if old_comp.num_atoms == struct.get("natoms"):
                composition: Composition = old_comp.copy()

            else: # The composition is a reduced form
                mult_factor: int = struct.get("natoms") // old_comp.num_atoms
                composition: Composition = old_comp * mult_factor

            composition: Dict[str, int] = composition.get_el_amt_dict()

        total_energy: float = struct.get("delta_e") * struct.get("natoms")
        new_name: str = "proc-oqmd-" + entry_id

        processed_data.update(
            {
                new_name: {
                    "entry_id": entry_id,
                    "composition": composition,
                    "final_energy": total_energy
                }
            }
        )
    new_filename = os.path.basename(path).replace(".json", "_proc.json")
    new_path = os.path.join(os.path.dirname(path), new_filename)

    with open(new_path, "wt") as fp:
        json.dump(processed_data, fp)


########################################


def load_phase_diagram_entries(filename: str) -> dict:
    """Load entries data from a JSON file."""
    print(f"Loading {filename}...")
    with open(filename, "r") as fp:
        entries = json.load(fp)

    for entry in entries.values():
        entry["composition"] = Composition.from_dict(entry["composition"])
    print(f"Loading finished, {len(entries)} entries found.")
    return entries


########################################

