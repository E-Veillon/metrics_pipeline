#!/usr/bin/python
"""Downloads and manages remote databases API and data."""

import os
import json
import typing as tp
import typing_extensions as tpe
from dataclasses import dataclass
import itertools as itt

from mp_api.client import MPRester
from emmet.core.types.enums import ThermoType
from emmet.core.thermo import ThermoDoc

from pymatgen.core import SETTINGS, Composition
from pymatgen.analysis.phase_diagram import PDEntry

from .io_base import PathLike
from .json import JsonLoader, JsonWriter

from metrics_pipeline.src.utils import ALL_ELT_SYMBOL_TO_Z

@dataclass
class PDEntryParser:
    """Parse structure data to get corresponding PDEntry object."""
    id: str
    composition: Composition | None = None
    formula: str | None = None
    natoms: int | None = None
    energy: float | None = None
    energy_per_atom: float | None = None

    def __post_init__(self) -> None:
        assert isinstance(self.id, str), TypeError(
            f"'id' expected a type 'str', got {type(self.id).__name__}."
        )

        # Get composition from formula and natoms
        if self.composition is None:
            assert self.formula is not None and self.natoms is not None, ValueError(
                f"Either 'composition' or 'formula' and 'natoms' must be given."
            )
            formula_unit = Composition(self.formula, strict=True)
            if formula_unit.num_atoms != self.natoms:
                mult_factor = self.natoms / formula_unit.num_atoms
                self.composition = formula_unit * mult_factor

        elif isinstance(self.composition, dict):
            self.composition = Composition(self.composition, strict=True)

        else:
            assert isinstance(self.composition, Composition), TypeError(
                "'composition' expected a type 'Composition', "
                f"got {type(self.composition).__name__}."
            )

        # Get total energy from energy per atom and composition atoms
        if self.energy is None:
            assert self.energy_per_atom is not None, ValueError(
                f"Either 'energy' or 'energy_per_atom' must be given."
            )
            assert self.composition is not None, "Type checker assertion."
            self.energy = self.energy_per_atom * self.composition.num_atoms

    @classmethod
    def from_composition_dict(cls, id: str, composition: dict[str, float], **kwargs) -> tpe.Self:
        """Build PDEntryParser with a composition dict instead of Composition object."""
        return cls(id, Composition(composition), **kwargs)

    def get_entry(self) -> PDEntry:
        """Build PDEntry from data."""
        assert isinstance(self.composition, Composition), "Type checker assertion."
        assert self.energy is not None, "Type checker assertion."

        return PDEntry(self.composition, self.energy, self.id)


class PDDataset:
    """Process structure database for phase diagram computations compatibility."""
    def __init__(
        self,
        data: dict,
        id_key: str = "entry_id",
        composition_key: str | None = None,
        formula_key: str | None = None,
        natoms_key: str | None = None,
        energy_key: str | None = None,
        energy_per_atom_key: str | None = None,
        attribute: str | None = None,
        compact: bool = False
    ) -> None:
        """
        Process structure database for phase diagram computations compatibility.

        Parameters
        ----------
        data: dict
            Database dict containing dicts of structures data to process.

        id_key: str, optional
            Mandatory key used to identify uniquely the structures in the database.
            Defaults to the key already compatible with phase diagrams, i.e. 'entry_id'.
            WARNING: If IDs are not unique, structures will override each other !

        composition_key: str, optional
            Key used to store structure composition as a dict of {"element": amount}.
            If not given, `formula_key` and `natoms_key` must be given.

        formula_key: str, optional
            Key used to store structure formula as a string. If not given,
            `composition_key` must be given.

        natoms_key: str, optional
            Key used to store the total number of atoms in unit cell. If not given,
            `composition_key` must be given.

        energy_key: str, optional
            Key used to store structure total energy in eV. If not given, `energy_per_atom_key`
            must be given.

        energy_per_atom: str, optional
            Key used to store structure energy per atom in eV/atom. If not given,
            `energy_key` must be given.

        attribute: str, optional
            Optional label to put on all entries in this dataset.

        compact: bool
            Whether loaded dataset is organized by data type instead of by structure.
            - If `False`, dataset is assumed to be organized by structure, i.e.
            {"s1": {s1_data}, "s2": {s2_data}, ...}.
            - If `True`, dataset is assumed to be organized by data type, i.e.
            {"entry_id": [s1, s2, ...], "composition": [s1, s2, ...], ...}.

            Defaults to `False`.
        """
        assert isinstance(data, dict), TypeError(
            f"'data' expected a type 'dict', got {type(data).__name__}."
        )
        val_types = ", ".join(sorted(set(type(val).__name__ for val in data.values())))
        assert all(isinstance(val, dict) for val in data.values()), TypeError(
            f"'data' values must contain only dict, got following types: {val_types}."
        )
        assert all(isinstance(val.get(id_key), str) for val in data.values()), TypeError(
            f"Given 'id_key' ({id_key}) does not always match a string ID in the data."
        )
        assert composition_key or (formula_key and natoms_key), ValueError(
            f"Either 'composition_key' or 'formula_key' and 'natoms_key' must be given."
        )
        assert energy_key or energy_per_atom_key, ValueError(
            f"Either 'energy_key' or 'energy_per_atom_key' must be given."
        )

        self._data: dict[str, PDEntry] = {}
        self._elements: set[str] = set()
        self._computable_elements: set[str] = set()
        self.key_dict = {
            "id": id_key,
            "composition": composition_key,
            "formula": formula_key,
            "natoms": natoms_key,
            "energy": energy_key,
            "energy_per_atom": energy_per_atom_key,
        }

        if compact:
            # Regenerate parsed structures data dicts
            parsed_data = {k: data.get(v, itt.repeat(None)) for k, v in self.key_dict.items()}
            keys = sorted(parsed_data.keys())
            data_list = [dict(zip(keys, values)) for values in zip(*(parsed_data[k] for k in keys))]

            # Iterate over already parsed structure data
            for struct_data in data_list:
                entry_parser = PDEntryParser(**struct_data)
                entry = entry_parser.get_entry()
                self._data[entry_parser.id] = entry
                self._elements.update(
                    set(entry.composition.get_el_amt_dict().keys())
                )
                if len(entry.composition) < 11:
                    self._computable_elements.update(
                        set(entry.composition.get_el_amt_dict().keys())
                    )

        else:
            # parse each structure data dict
            for _, struct_data in data.items():
                processed_data = {k: struct_data.get(v) for k, v in self.key_dict.items()}
                entry_parser = PDEntryParser(**processed_data) # type: ignore
                entry = entry_parser.get_entry()
                self._data[entry_parser.id] = entry
                self._elements.update(
                    set(entry.composition.get_el_amt_dict().keys())
                )
                if len(entry.composition) < 11:
                    self._computable_elements.update(
                        set(entry.composition.get_el_amt_dict().keys())
                    )

        if attribute is not None:
            for entry in self._data.values():
                entry.attribute = attribute

    def get_all_entries(self) -> dict[str, PDEntry]:
        """Get dataset dict of all phase diagram entries."""
        return self._data

    def get_computable_entries(self) -> dict[str, PDEntry]:
        """
        Get dataset dict of entries with element dimension at most 10, which is the maximum
        supported by the phase diagram constructor.
        """
        return {name: entry for name, entry in self._data.items() if len(entry.composition) < 11}

    def get_uncomputable_entries(self) -> dict[str, PDEntry]:
        """
        Get dataset dict of entries with element dimension greater than 10, which is the maximum
        supported by the phase diagram constructor.
        """
        return {name: entry for name, entry in self._data.items() if len(entry.composition) > 10}

    def get_filtered_entries(
        self, elts: set[str] | None = None, dims: set[int] | None = None
    ) -> dict[str, PDEntry]:
        """
        Flexible query method to get filtered dataset dict of entries.

        Parameters
        ----------
        elts: list[str], optional
            List of elements the entries compositions must fit in.
            Only entries containing only elements in the list are returned.
            If not given, no restriction is applied.

        dims: list[int], optional
            List of accepted element dimensions.
            Only entries having a composition of these dimensions will be returned.
            If not given, no restriction is applied.

        Returns
        -------
        dict[str, PDEntry]
            Dataset dict of entries corresponding to given conditions.
        """
        if elts is None and dims is None:
            return self.get_all_entries()

        data_tuples = list(self._data.items())

        if elts is not None:
            data_tuples = list(
                filter(
                    lambda data_tup: all(elt in elts for elt in data_tup[1].elements),
                    data_tuples
                )
            )
        if dims is not None:
            data_tuples = list(
                filter(
                    lambda data_tup: len(data_tup[1].composition) in dims,
                    data_tuples
                )
            )
        return dict(data_tuples)

    @property
    def max_dim(self) -> int:
        """Max number of distinct elements in all entries."""
        return max(len(entry.elements) for entry in self.get_all_entries().values())

    @property
    def max_computable_dim(self) -> int:
        """Max number of distinct elements in computable entries."""
        return max(len(entry.elements) for entry in self.get_computable_entries().values())

    @property
    def elements(self) -> set[str]:
        """Set of elements symbols for all elements used in the dataset."""
        return self._elements

    @property
    def elements_alphabetic(self) -> list[str]:
        """
        List of unique elements symbols for all elements used in the dataset.
        Sorted alphabetically.
        """
        return sorted(self._elements)

    @property
    def elements_periodic(self) -> list[str]:
        """
        List of unique elements symbols for all elements used in the dataset.
        Sorted by atomic number.
        """
        return sorted(self._elements, key=lambda symbol: ALL_ELT_SYMBOL_TO_Z[symbol])

    @property
    def computable_elements(self) -> set[str]:
        """Set of element symbols for elements used in computable entries."""
        return self._computable_elements

    @property
    def computable_elements_alphabetic(self) -> list[str]:
        """
        List of unique elements symbols for elements used in computable entries.
        Sorted alphabetically.
        """
        return sorted(self._computable_elements)

    @property
    def computable_elements_periodic(self) -> list[str]:
        """
        List of unique elements symbols for elements used in computable entries.
        Sorted by atomic number.
        """
        return sorted(self._computable_elements, key=lambda symbol: ALL_ELT_SYMBOL_TO_Z[symbol])

    @classmethod
    def from_file(cls, filepath: PathLike, **kwargs) -> tpe.Self:
        """
        Load the dataset from a JSON file instead of initializing with a dict.

        Parameters
        ----------
        filepath: str | Path
            Path to the JSON dataset file.

        kwargs: Any
            Any keyword argument to pass to the constructor.
        """
        loaded_data = JsonLoader(filepath).load_as_dict()
        return cls(loaded_data, **kwargs)

    def write_dataset(self, filepath: PathLike, compact: bool = False) -> None:
        """
        Write processed dataset to a JSON file.

        Parameters
        ----------
        filepath: str | Path
            Path to the JSON file to write.

        compact: bool
            Whether to reduce overhead by organizing by data type instead of by structure.

            - If `False`, data is organized by structure, i.e.
            {"s1": {s1_data}, "s2": {s2_data}, ...}.
            - If `True`, data is organized by data type, i.e.
            {"entry_id": [s1, s2, ...], "composition": [s1, s2, ...], ...}.

            This method is more memory efficient as soon as the file contains more structures than
            data types. Defaults to False.
        """
        if compact:
            # Ensure data ordering with a list
            data_list = list(self._data.items())
            data = {
                "entry_id": [entry_id for entry_id, _ in data_list],
                "composition": [entry.composition.get_el_amt_dict() for _, entry in data_list],
                "final_energy": [entry.energy for _, entry in data_list]
            }
        else:
            data = {
                entry_id: {
                    "entry_id": entry_id,
                    "composition": entry.composition.get_el_amt_dict(),
                    "final_energy": entry.energy
                } for entry_id, entry in self._data.items()
            }
        JsonWriter(filepath, data).write_as_dict()


class APIKeyNotFoundError(Exception):
    """
    An Error ocurring when an API key is not defined while trying to access an API.
    """


class MPDatasetDownloader:
    """
    Download the Materials Project phase diagram dataset and format the data
    for compatibility with local phase diagram computations.
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
        # Used for file assertions
        self.writer = JsonWriter(save_file)
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
            data = tp.cast(list[dict], data)
            dataset = PDDataset(
                {entry["material_id"]: entry for entry in data},
                id_key="material_id",
                composition_key="composition",
                energy_per_atom_key="energy_per_atom"
            )
            dataset.write_dataset(self.writer.filepath)

    def print_avail_fields(self, endpoint: str) -> None:
        """Print the list of available fields corresponding to given MPRester endpoint."""
        with MPRester(self.api_key) as mpr:
            match endpoint:
                case "materials":
                    print(mpr.materials.available_fields)
                case "thermo":
                    print(mpr.materials.thermo.available_fields)
                case _:
                    raise NotImplementedError(
                        f"Given 'endpoint' ({endpoint}) does not match any implemented doc."
                    )
