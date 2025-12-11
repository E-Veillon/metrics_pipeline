#!/usr/bin/python
"""Downloads and manages remote databases API and data."""

import os
import json
from typing import Dict, Any
from mp_api.client import MPRester
from emmet.core.types.enums import ThermoType
from emmet.core.thermo import ThermoDoc

from pymatgen.core import SETTINGS, Composition

from .io_base import check_file_format, check_file_or_dir, PathLike


class APIKeyNotFoundError(Exception):
    """
    An Error ocurring when an API key is not defined while trying to access an API.
    """


class MPDatasetDownloader:
    """
    Download Materials Project dataset and format the data for compatibility with GenMat.
    """
    def __init__(
        self, save_file: PathLike, api_key: str | None = None, download: bool = True
    ) -> None:
        """
        Download Materials Project dataset and format the data for compatibility with GenMat.

        Parameters
        ----------
        save_file: str | Path
            Path to a file to save the formatted dataset. Must have JSON file format.

        api_key: str, optional
            Give manually your MP API key. If not given, will retrieve it from `PMG_MAPI_KEY`
            in `.pmgrc.yaml`.

        download: bool
            Whether to download the dataset at initialization. Set it to False if you just need
            to get MP API infos needing an API key without actual download. Defaults to True.
        """
        check_file_format(save_file, allowed_formats="json")
        os.makedirs(os.path.dirname(save_file), exist_ok=True)
        self.save_file = save_file
        self.api_key = self._get_api_key(api_key)

        if download:
            self.download()

    @staticmethod
    def _get_api_key(api_key: str | None = None) -> str:
        """
        Check whether the MP API key is defined either in .pmgrc.yaml or manually.
        If it is found, the key string is returned. Otherwise, an exception is raised.
        """
        if api_key is not None:
            assert isinstance(api_key, str), TypeError(
                f"'api_key' expected a type 'str', got {type(api_key).__name__}."
            )
        api_key = api_key or SETTINGS.get("PMG_MAPI_KEY")
        if not api_key:
            raise APIKeyNotFoundError(
                "An API key must be given to access MP API. "
                "Configure it in .pmgrc.yaml in 'PMG_MAPI_KEY' keyword or pass it manually."
            )
        return api_key

    def download(self) -> None:
        """
        Download GGA/GGA+U phase diagram relevant entries from MP API
        to a local JSON file. Downloaded entries are dicts of the form
        {"entry_id": str, "composition": dict, "final_energy": float}.
        """
        with MPRester(
            api_key=self.api_key,
            use_document_model=False
            ) as mpr:
            data = mpr.materials.thermo.search(thermo_types=[ThermoType.GGA_GGA_U],
                all_fields=False, fields=["material_id","composition","energy_per_atom"]
            )
            print(f"Number of Materials Project entries downloaded: {len(data)}")

            final_data = {}
            for entry in data:
                assert isinstance(entry, dict), RuntimeError(
                    "Should never trigger, but satisfies type checker."
                )
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

        with open(self.save_file, "wt", encoding="utf-8") as fp:
            json.dump(final_data, fp)

    def print_avail_fields(self, endpoint: str) -> None:
        """Print the list of available fields corresponding to given MPRester endpoint."""
        with MPRester(self.api_key) as mpr:
            if endpoint == "materials":
                print(mpr.materials.available_fields)
            elif endpoint == "thermo":
                print(mpr.materials.thermo.available_fields)
            else:
                raise NotImplementedError(
                    "Given 'endpoint' arg does not match any implemented doc."
                )


def process_oqmd_json_file(path: str) -> None:
    """
    Process OQMD json data file to get needed data format in a new JSON file.
    Expected OQMD JSON key formatting is the one got with OQMD RESTful API.
    The processed file have a '_proc' suffix added to its name.

    Parameters:
        path (str): Path to the OQMD JSON file to process.
    """
    check_file_or_dir(path, "file", allowed_formats="json")

    with open(path, "rt", encoding="utf-8") as fp:
        data = json.load(fp)

    processed_prefix = "proc-oqmd-"
    processed_data = {}

    for name, struct in data.items():
        struct: Dict[str, Any]

        if (
            struct.get("entry_id") is None
            or struct.get("composition") is None
            or struct.get("natoms") is None
            or struct.get("delta_e") is None
        ):
            if name.startswith(processed_prefix):
                raise ValueError(
                    f"structure '{name}' was already processed (begins with '{processed_prefix}'). "
                    "Please verify that you are not processing the same data more than "
                    "once, as it may cause wrong total energy computations."
                )
            raise ValueError(
                f"At least one of needed data keys is missing in structure '{name}'. "
                "Please verify that the following keys are present: "
                "'entry_id', 'composition', 'natoms', 'delta_e'."
            )

        entry_id: str = str(struct.get("entry_id"))

        if isinstance(struct.get("composition"), Dict):
            comp_dict: Dict[str, int|float] = struct["composition"]

        elif isinstance(struct["composition"], str):
            old_comp = Composition("".join(struct["composition"].split()), strict=True)

            if old_comp.num_atoms == struct["natoms"]:
                composition: Composition = old_comp.copy()

            else: # The composition is a reduced form
                mult_factor: int = struct["natoms"] // old_comp.num_atoms
                composition: Composition = old_comp * mult_factor

            comp_dict: Dict[str, int|float] = composition.get_el_amt_dict()

        total_energy: float = struct["delta_e"] * struct.get("natoms")
        new_name: str = processed_prefix + entry_id

        processed_data.update(
            {
                new_name: {
                    "entry_id": entry_id,
                    "composition": comp_dict,
                    "final_energy": total_energy
                }
            }
        )
    new_filename = os.path.basename(path).replace(".json", "_proc.json")
    new_path = os.path.join(os.path.dirname(path), new_filename)

    with open(new_path, "wt", encoding="utf-8") as fp:
        json.dump(processed_data, fp)


def load_phase_diagram_entries(filename: str) -> dict:
    """Load entries data from a JSON file."""
    check_file_or_dir(filename, "file", allowed_formats="json")

    print(f"Loading {filename}...")
    with open(filename, "rt", encoding="utf-8") as fp:
        entries = json.load(fp)

    for entry in entries.values():
        entry["composition"] = Composition.from_dict(entry["composition"])
    print(f"Loading finished, {len(entries)} entries found.")
    return entries
