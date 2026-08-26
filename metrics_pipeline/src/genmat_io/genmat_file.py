"""
Class handling specific JSON files containing GenMat structures data in a predictable and compact
way to avoid ambiguities in data communications between pipeline steps.
"""

from __future__ import annotations
import typing as tp
import typing_extensions as tpe
from collections.abc import Callable
import itertools as itt

from pymatgen.core import Structure

from .structure_file import (
    StructureFile, _parse_dict_keys,
    _compress_struct_dict, _uncompress_struct_dict
)
from src.utils.genmat_data import GenMatStructure
from src.utils.spg_data import Spacegroup


G = tp.TypeVar("G", bound=GenMatStructure)


class GenMatFile(StructureFile[G]):
    """
    I/O class to load and write JSON files containing structure data for the GenMat library.
    """
    # Extension keys
    _genmat_mandatory_keys: tuple[str, ...] = ("name",)
    _genmat_optional_keys: tuple[str, ...] = (
        "special_keys", "spacegroup", "energy", "energy_above_hull"
    )
    _genmat_keys: tuple[str, ...] = _genmat_mandatory_keys + _genmat_optional_keys

    # Extension of dict representations keys
    _dict_mandatory_keys: tuple[str, ...] = (
        _genmat_mandatory_keys + StructureFile._dict_mandatory_keys
    )
    _dict_optional_keys: tuple[str, ...] = (
        _genmat_optional_keys + StructureFile._dict_optional_keys
    )

    # Extension of internal storage-efficient keys
    _file_mandatory_keys: tuple[str, ...] = (
        _genmat_mandatory_keys + StructureFile._file_mandatory_keys
    )
    _file_optional_keys: dict[str, Callable[[tpe.Self], tp.Iterator[tp.Any]]] = {
        "special_keys": lambda self: iter({} for _ in range(self.length)),
        "spacegroup": lambda self: iter(Spacegroup(0) for _ in range(self.length)),
        "energy": lambda self: itt.repeat(None, self.length),
        "energy_above_hull": lambda self: itt.repeat(None, self.length)
    } | StructureFile._file_optional_keys

    def _build_structure(self, index: int) -> G:
        """Use data at given index to build corresponding GenMatStructure."""
        dct = _uncompress_struct_dict(
            {k: self._data[k][index] for k in self.file_active_keys}
        )
        return GenMatStructure.from_dict(dct)

    def parse_pmg_structures(self) -> list[Structure]:
        """Build all standard pymatgen `Structure` objects from file data, in order."""
        return [gstruct.as_structure() for gstruct in self]

    @property
    def all_names_iter(self) -> tp.Iterator[str]:
        """Lazy iterator of all stored structure names."""
        return iter(self._data["name"])

    @property
    def all_special_keys_iter(self) -> tp.Iterator[dict[str, tp.Any]]:
        """Lazy iterator of all stored special keys dictionaries."""
        if "special_keys" in self.saved_optional_keys:
            return iter(self._data["special_keys"])
        return self._get_optional_key_default("special_keys")

    @property
    def all_spacegroups(self) -> tp.Iterator[Spacegroup]:
        """
        Lazy iterator of all stored spacegroups, converted into Spacegroup objects on the fly.
        """
        if "spacegroup" in self.saved_optional_keys:
            return iter(Spacegroup.from_dict(dct) for dct in self._data["spacegroup"])
        return self._get_optional_key_default("spacegroup")

    @property
    def all_energies_iter(self) -> tp.Iterator[float | None]:
        """Lazy iterator of all stored energy values."""
        if "energy" in self.saved_optional_keys:
            return iter(self._data["energy"])
        return self._get_optional_key_default("energy")

    @property
    def all_energies_above_hull_iter(self) -> tp.Iterator[float | None]:
        """Lazy iterator of all stored energy above hull values."""
        if "energy_above_hull" in self.saved_optional_keys:
            return iter(self._data["energy_above_hull"])
        return self._get_optional_key_default("energy_above_hull")
