"""
Class handling specific JSON files containing structures data in a predictable and compact
way to avoid ambiguities in data communications between pipeline steps.
"""

import typing as tp
import typing_extensions as tpe
from collections.abc import Mapping, Callable
import functools as ft
import itertools as itt
from enum import IntEnum

from tqdm.contrib.concurrent import process_map
import numpy as np

from pymatgen.core import Structure, Lattice, PeriodicSite

from .io_base import PathLike
from .json import JsonLoader, JsonWriter
from src.utils.common_asserts import check_type
from src.utils.visual_iterator import VisualIterator


class FloatPrecision(IntEnum):
    """Enum of possible float values precisions."""
    LOW = 4
    MEDIUM = 8
    HIGH = 16


class StructureFile:
    """
    I/O class to load and write JSON files containing efficiently stored structure data.
    """
    _file_keys: tuple[str, ...] = (
        "a", "b", "c", "alpha", "beta", "gamma", "charge", "properties",
        "species", "site_a", "site_b", "site_c", "site_properties"
    )
    _lattice_keys: tuple[str, ...] = (
        "a", "b", "c", "alpha", "beta", "gamma"
    )
    _site_keys: tuple[str, ...] = (
        "species", "site_a", "site_b", "site_c", "site_properties"
    )
    _dict_keys: tuple[str, ...] = (
        "lattice", "charge", "properties", "sites"
    )
    _optional_keys: dict[str, Callable[[tpe.Self], tp.Iterator[tp.Any]]] = {
        "charge": lambda self: itt.repeat(0.0, self.length), # Fastest, for immutable defaults
        "properties": lambda self: iter({} for _ in range(self.length)), # mutable defaults
        "site_properties": lambda self: iter( # site-wise attribute defaults
            [{} for _ in range(n_sites)] for n_sites in self.num_sites_per_struct_iter
        )
    }
    _data: dict[str, list[tp.Any]]
    _float_prec: FloatPrecision

    def __init__(
        self,
        data: dict[str, list[tp.Any]] | None = None,
        compressed: bool = False,
        float_prec: FloatPrecision | str | int = FloatPrecision.HIGH,
        saved_optional_keys: set[str] | None = None
    ) -> None:
        """
        I/O class to load and write JSON files containing efficiently stored structure data.

        Parameters
        ----------
        data: dict[str, list[Any]], optional
            Structure data to put into the file. Must be a dict with the same keys
            as a dict generated with `Structure.as_dict()`, containing lists of equal
            lengths where the same index contains data for the same structure. If not given,
            instantiate an empty file that can be filled in via `add_structure()` and
            `add_structures()` methods.

        compressed: bool
            Whether given data already follows the storage-efficient dict formatting, which
            optimizes the `structure.as_dict()` keys into several easier-to-store values.
            Defaults to False.

        float_prec: FloatPrecision | str | int
            What level of precision to apply to lattice parameterss and site positions
            values when writing structure data into a file:
            - 'low' or 4: rounds to 4 decimal places
            (0.01 pm for lengths and positions, 1e-4 degree for angles);
            - 'medium' or 8: rounds to 8 decimal places;
            - 'high' or 16: rounds to 16 decimal places (numerical precision level).

            WARNING: Lost precision is irreversible and will NOT be recovered when retrieving data.
            Defaults to 'high' to avoid unintended precision loss.

            Can be modified after initialization, and be applied explicitly on stored data without
            writing to file by using the `.round()` method.

        saved_optional_keys: set[str], optional
            Pass optional data attributes that should be stored in this file. If not given,
            all optional attributes are stored by default. Use this to shave off attributes that
            are not used or computed in this file, increasing storage efficiency by avoiding
            storage of empty attributes.
        """
        self.float_prec = float_prec

        if saved_optional_keys is None:
            saved_optional_keys = set(self._optional_keys)
        else:
            for key in saved_optional_keys:
                self._get_optional_key_default(key) # test optional key existence

        self._saved_optional_keys = saved_optional_keys

        if data is None:
            self._data = {}
            for key in self.active_file_keys:
                self._data[key] = []
            return

        keys_to_check = self.active_file_keys if compressed else self._dict_keys
        _check_dict_keys(data, keys_to_check)

        if not len(len_set:={len(data[key]) for key in keys_to_check}) == 1:
            len_set_str = ", ".join(sorted(map(str, len_set)))
            raise ValueError(
                f"Some iterables in data do not have the same length, got {len_set_str}."
            )
        
        if not compressed:
            # Reorganize keys for storage-efficient format
            ordered_data = (data[k] for k in self._dict_keys)
            parsed_data = (
                _compress_struct_dict(
                    {k: v for k, v in zip(self._dict_keys, struct_data, strict=True)}
                ) for struct_data in zip(*ordered_data, strict=True)
            )
            compressed_data = {key: [] for key in self.active_file_keys}
            for dct in parsed_data:
                for key in self.active_file_keys:
                    compressed_data[key].append(dct[key])

            data = compressed_data

        # Handle saved/unsaved optional keys
        for opt_key in self._optional_keys:
            if opt_key in self.saved_optional_keys and opt_key not in data:
                data[opt_key] = list(self._get_optional_key_default(key))

            elif opt_key not in self.saved_optional_keys:
                data.pop(key, None)

        self._data = data

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, index: int | slice) -> Structure | tpe.Self:
        """Get a Structure from the data at given index, or a sub-dataset from a slice."""
        if isinstance(index, slice):
            subdata = {key: self._data[key][index] for key in self.active_file_keys}
            cls = type(self)
            return cls(
                subdata,
                compressed=True,
                saved_optional_keys=self.saved_optional_keys
            )
        elif isinstance(index, int):
            if index < 0:
                index += self.length
            if not 0 <= index < self.length:
                raise IndexError(f"{type(self).__name__}: index out of range.")

            return self._build_structure(index)

        check_type(index, "index", (int, slice))

    def __iter__(self) -> tp.Iterator[Structure]:
        return iter(self._build_structure(i) for i in range(self.length))

    def __add__(self, other: tpe.Self) -> tpe.Self:
        cls = type(self)

        if not isinstance(other, cls):
            return NotImplemented

        new_dataset = {}

        # Handle data concatenation with possibly lacking optional keys
        all_active_file_keys = set(self.active_file_keys) | set(other.active_file_keys)
        self_lacking_keys = other.saved_optional_keys - self.saved_optional_keys
        other_lacking_keys = self.saved_optional_keys - other.saved_optional_keys

        for key in all_active_file_keys: # No iteration over optional keys unsaved on both sides
            if key in self_lacking_keys: # Optional key not in self but present in other
                new_dataset[key] = list(self._get_optional_key_default(key)) + other._data[key]
            elif key in other_lacking_keys: # Optional key present in self but not in other
                new_dataset[key] = self._data[key] + list(other._get_optional_key_default(key))
            else: # Mandatory or optional key present in both
                new_dataset[key] = self._data[key] + other._data[key]

        # Keep highest float precision to avoid unintended truncature
        new_float_prec = max(self.float_prec, other.float_prec)

        # Merge active optional keys from both datasets
        new_saved_optional_keys = self.saved_optional_keys | other.saved_optional_keys

        return cls(
            new_dataset,
            compressed=True,
            float_prec=new_float_prec,
            saved_optional_keys=new_saved_optional_keys
        )

    def _build_structure(self, index: int) -> Structure:
        """Use data at given index to build corresponding Structure."""
        dct = _uncompress_struct_dict(
            {k: self._data[k][index] for k in self.active_file_keys}
        )
        return Structure.from_dict(dct)

    def _check_ordering(self, structure: Structure) -> None:
        """Check that passed structure is ordered."""
        if not structure.is_ordered:
            raise ValueError(
                f"{type(self).__name__!r} only supports ordered structures, "
                "i.e. with only one species with occupancy 1.0 per site."
            )

    def _get_optional_key_default(self, key: str) -> tp.Iterator[tp.Any]:
        """Get an iterator containing default values for an optional data key."""
        if key not in self._optional_keys:
            raise ValueError(f"{key!r} is not a valid optional key.")
        return self._optional_keys[key](self)

    def copy(self) -> tpe.Self:
        """Get a shallow copy of this file instance."""
        cls = type(self)
        return cls(
            data=self._data,
            compressed=True,
            float_prec=self.float_prec,
            saved_optional_keys=self.saved_optional_keys
        )

    def round(self) -> tpe.Self:
        """
        Round all stored lattice parameters and site positions using stored `float_prec`.
        This operation is irreversible and is done in place.
        """
        # Use numpy vectorized rounding for lattice parameters lists
        lattice_array = np.array(
            [self._data[key] for key in self._lattice_keys], dtype=np.float64
        ).T
        np.round(lattice_array, decimals=self.float_prec.value, out=lattice_array)
        rounded_lists = lattice_array.T.tolist()

        for key, rounded_list in zip(self._lattice_keys, rounded_lists):
            self._data[key] = rounded_list

        # Use python loop for site-wise rounding
        site_coords_iter = zip(
            self.all_site_a_coords_iter,
            self.all_site_b_coords_iter,
            self.all_site_c_coords_iter
        )
        for idx, (site_a_list, site_b_list, site_c_list) in enumerate(site_coords_iter):
            self._data["site_a"][idx] = [
                round(a, self.float_prec) for a in site_a_list
            ]
            self._data["site_b"][idx] = [
                round(b, self.float_prec) for b in site_b_list
            ]
            self._data["site_c"][idx] = [
                round(c, self.float_prec) for c in site_c_list
            ]

        return self

    def parse_structures(self, verbose: bool = False) -> list[Structure]:
        """
        Build all Structure objects from file data, in order.

        Parameters
        ----------
        verbose: bool
            Whether to show parsing progression with a VisualIterator.
            Defaults to False.

        Returns
        -------
        list[Structure]
            List of Structure objects extracted from data.
        """
        iterator = iter(self)

        if verbose:
            desc = "Extracting structures from file"
            iterator = VisualIterator.from_big_iterator(
                iterator, n_elts=self.length, desc=desc, unit="extracted", percent=True
            )

        return list(iterator)

    def add_structure(self, structure: Structure) -> None:
        """Add one Structure object data to the file."""
        self._check_ordering(structure)
        dct = structure.as_dict(verbosity=0)
        dct = _compress_struct_dict(dct)
        _check_dict_keys(dct, self.active_file_keys)
        for key in self.active_file_keys:
            self._data[key].append(dct[key])

    def add_structures(self, structures: tp.Iterable[Structure]) -> None:
        """
        Convenient method to add a full iterable of structure data to the file.
        Order is preserved.
        """
        for structure in structures:
            self.add_structure(structure)

    @property
    def saved_optional_keys(self) -> set[str]:
        """Set of active optional keys."""
        return self._saved_optional_keys

    @property
    def active_file_keys(self) -> tuple[str, ...]:
        """Get internal keys where data is actively saved."""
        file_keys = list(self._file_keys)

        for opt_key in self._optional_keys:
            if opt_key not in self.saved_optional_keys:
                file_keys.remove(opt_key)

        return tuple(file_keys)

    @property
    def length(self) -> int:
        """Number of structures contained in the file."""
        assert len(len_set:={len(self._data[key]) for key in self.active_file_keys}) == 1, (
            f"{type(self).length.__qualname__}: The lengths of data lists in the file "
            "are not equal, some structure data must be incomplete or corrupted."
        )
        return len_set.pop()

    @property
    def all_a_params_iter(self) -> tp.Iterator[float]:
        """Lazy iterator of all stored lattice parameter 'a' in order."""
        return iter(self._data["a"])

    @property
    def all_b_params_iter(self) -> tp.Iterator[float]:
        """Lazy iterator of all stored lattice parameter 'b' in order."""
        return iter(self._data["b"])

    @property
    def all_c_params_iter(self) -> tp.Iterator[float]:
        """Lazy iterator of all stored lattice parameter 'c' in order."""
        return iter(self._data["c"])

    @property
    def all_alpha_angles_iter(self) -> tp.Iterator[float]:
        """Lazy iterator of all stored lattice parameter 'alpha' in order."""
        return iter(self._data["alpha"])

    @property
    def all_beta_angles_iter(self) -> tp.Iterator[float]:
        """Lazy iterator of all stored lattice parameter 'beta' in order."""
        return iter(self._data["beta"])

    @property
    def all_gamma_angles_iter(self) -> tp.Iterator[float]:
        """Lazy iterator of all stored lattice parameter 'gamma' in order."""
        return iter(self._data["gamma"])

    @property
    def all_lattices_iter(self) -> tp.Iterator[Lattice]:
        """Lazy iterator of all Lattice objects rebuilt from data, in order."""
        params_iter = zip(*(self._data[key] for key in self._lattice_keys))
        return (
            Lattice.from_parameters(a, b, c, alpha, beta, gamma)
            for a, b, c, alpha, beta, gamma in params_iter
        )

    @property
    def all_charges_iter(self) -> tp.Iterator[float]:
        """
        Lazy iterator of all stored structure charges in order. If `save_charge` was set
        to False, builds an iterator of default charges of 0.0 coherent with file data
        for file iteration consistency.
        """
        if "charge" in self.saved_optional_keys:
            return iter(self._data["charge"])
        return self._get_optional_key_default("charge")

    @property
    def all_properties_iter(self) -> tp.Iterator[dict]:
        """
        Lazy iterator of all stored structure properties dicts in order. If `save_properties`
        was set to False, builds an iterator of empty properties dicts coherent with file data
        for file iteration consistency.
        """
        if "properties" in self.saved_optional_keys:
            return iter(self._data["properties"])
        return self._get_optional_key_default("properties")

    @property
    def all_species_iter(self) -> tp.Iterator[list[str]]:
        """
        Lazy iterator of all stored lists of species (one list per structure)
        in order.
        """
        return iter(self._data["species"])

    @property
    def num_sites_per_struct_iter(self) -> tp.Iterator[int]:
        """Lazy iterator of the number of sites in each stored structure."""
        return (len(species) for species in self.all_species_iter)

    @property
    def all_site_a_coords_iter(self) -> tp.Iterator[list[float]]:
        """
        Lazy iterator of all stored lists of sites fractional coordinates
        in the direction 'a' (one list per structure) in order.
        """
        return iter(self._data["site_a"])

    @property
    def all_site_b_coords_iter(self) -> tp.Iterator[list[float]]:
        """
        Lazy iterator of all stored lists of sites fractional coordinates
        in the direction 'b' (one list per structure) in order.
        """
        return iter(self._data["site_b"])

    @property
    def all_site_c_coords_iter(self) -> tp.Iterator[list[float]]:
        """
        Lazy iterator of all stored lists of sites fractional coordinates
        in the direction 'c' (one list per structure) in order.
        """
        return iter(self._data["site_c"])

    @property
    def all_site_properties_iter(self) -> tp.Iterator[list[dict]]:
        """
        Lazy iterator of all stored lists of site properties dicts (one list per structure)
        in order. If `save_site_properties` was set to False, builds an iterator of lists
        of empty site properties dicts coherent with file data for file iteration consistency.
        """
        if "site_properties" in self.saved_optional_keys:
            return iter(self._data["site_properties"])
        return self._get_optional_key_default("site_properties")

    @property
    def all_sites_iter(self) -> tp.Iterator[list[PeriodicSite]]:
        """
        Lazy iterator of all lists of sites (one list per structure)
        rebuilt from data, in order.
        """
        structs_iter = zip(*(self._data[key] for key in self._site_keys))
        lattices_iter = self.all_lattices_iter

        for sites_data_tuple, lattice in zip(structs_iter, lattices_iter):
            sites_iter = zip(*sites_data_tuple)
            yield [
                PeriodicSite(
                    elt, [site_a, site_b, site_c], lattice, site_props
                ) for elt, site_a, site_b, site_c, site_props in sites_iter
            ]

    @property
    def float_prec(self) -> FloatPrecision:
        """
        Initialized float precision level. Only applied when writing data to file
        or when the `.round()` method is used.
        """
        return self._float_prec

    @float_prec.setter
    def float_prec(self, value: FloatPrecision | str | int) -> None:
        if isinstance(value, FloatPrecision):
            self._float_prec = value
        elif isinstance(value, str) and value.upper() in {"LOW", "MEDIUM", "HIGH"}:
            self._float_prec = FloatPrecision[value.upper()]
        elif isinstance(value, int) and value in {4, 8, 16}:
            self._float_prec = FloatPrecision(value)
        else:
            raise ValueError(
                f"{value!r} is not a valid value for 'float_prec'. "
                f"Valid values are: a {FloatPrecision.__name__} object, "
                "the strings 'low', 'medium', or 'high' (case insensitive), "
                "and the integers 4, 8, or 16."
            )

    @classmethod
    def from_file(
        cls, filepath: PathLike, match_optional_keys: bool = True, **kwargs: tp.Any
    ) -> tpe.Self:
        """
        Load a JSON file of structure data.

        Parameters
        ----------
        filepath: str | Path
            Path to the file to load.

        match_optional_keys: bool
            Whether to build "saved_optional_keys" arg automatically to match keys that are
            present in file. Defaults to True. Avoids generating default empty data
            at initialization by default whenever loading a file. Set to False to generate
            default values for all unsaved optional keys. Manually passed "saved_optional_keys"
            keyword argument will override this argument behavior.

        kwargs: Any
            Additional keyword arguments to pass to the constructor, except "compressed"
            (auto set to True).
        """
        file_data = JsonLoader(filepath).load_as_dict()
        if file_data.get("@module") != cls.__module__:
            raise ValueError(
                f"Invalid module in file: expected {cls.__module__!r}, "
                f"got {file_data.get('@module')!r}"
            )
        if file_data.get("@class") != cls.__name__:
            raise ValueError(
                f"Invalid class in file: expected {cls.__name__!r}, "
                f"got {file_data.get('@class')!r}"
            )
        check_type(data:=file_data.get("data"), "data", (dict,))
        data = tp.cast(dict, data)

        kwargs["compressed"] = True

        if match_optional_keys:
            auto_saved_keys = set(filter(lambda key: key in data, cls._optional_keys))
            kwargs.setdefault("saved_optional_keys", auto_saved_keys)

        return cls(data, **kwargs)

    def write_file(self, filepath: PathLike) -> None:
        """Write stored data to a compact JSON file."""
        written_data = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "data": self.round()._data
        }
        JsonWriter(filepath, written_data).write_as_dict()


def _check_dict_keys(dct: dict[str, tp.Any], keys_to_check: tp.Sequence[str]) -> None:
    """Verify that all given keys are present in the dict."""
    if not all(keys_in_data:=[key in dct for key in keys_to_check]):
        lacking_keys = ", ".join([
            keys_to_check[idx]
            for idx, key_in_data in enumerate(keys_in_data) if not key_in_data
        ])
        raise KeyError(f"Some mandatory keys are not present in 'data': {lacking_keys}.")


def _compress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into more storage-efficient organization.
    Exact reverse operation of `_uncompress_struct_dict()`.
    """
    file_dict = {}
    file_dict["charge"] = dct.pop("charge")
    file_dict["properties"] = dct.pop("properties")

    # Flatten lattice parameters
    lattice = Lattice.from_dict(dct["lattice"])
    file_dict["a"] = lattice.a
    file_dict["b"] = lattice.b
    file_dict["c"] = lattice.c
    file_dict["alpha"] = lattice.alpha
    file_dict["beta"] = lattice.beta
    file_dict["gamma"] = lattice.gamma

    # Flatten sites data
    file_dict["species"] = []
    file_dict["site_a"], file_dict["site_b"], file_dict["site_c"] = [], [], []
    file_dict["site_properties"] = []
    for site_dict in dct["sites"]:
        file_dict["species"].append(site_dict["species"][0]["element"])
        file_dict["site_a"].append(site_dict["abc"][0])
        file_dict["site_b"].append(site_dict["abc"][1])
        file_dict["site_c"].append(site_dict["abc"][2])
        file_dict["site_properties"].append(site_dict["properties"])

    return file_dict


def _uncompress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into object reconstruction-friendly organization.
    Exact reverse operation of `_compress_struct_dict()`.
    """
    struct_dict = {}
    struct_dict["charge"] = dct.pop("charge", 0.0)
    struct_dict["properties"] = dct.pop("properties", {})

    # Rebuild pymatgen-style lattice dict
    struct_dict["lattice"] = Lattice.from_parameters(
        a=dct.pop("a"), b=dct.pop("b"), c=dct.pop("c"),
        alpha=dct.pop("alpha"), beta=dct.pop("beta"), gamma=dct.pop("gamma")
    ).as_dict()

    # Rebuild pymatgen-style site dict
    species_list = dct.pop("species")
    site_iter = zip(
        species_list,
        dct.pop("site_a"),
        dct.pop("site_b"),
        dct.pop("site_c"),
        dct.pop("site_properties", [{} for _ in range(len(species_list))]),
        strict=True
    )
    struct_dict["sites"] = [
        {
            "species": [{"element": site_data[0], "occu": 1.0}],
            "abc": [site_data[1], site_data[2], site_data[3]],
            "properties": site_data[4]
        } for site_data in site_iter
    ]
    return struct_dict
