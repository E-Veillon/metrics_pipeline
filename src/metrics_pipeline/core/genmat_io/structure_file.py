"""
Class handling specific JSON files containing structures data in a predictable and compact
way to avoid ambiguities in data communications between pipeline steps.
"""

from __future__ import annotations
import typing as tp
import typing_extensions as tpe
from collections.abc import Callable
import itertools as itt

import numpy as np

from pymatgen.core import Structure, Lattice, PeriodicSite

from .io_base import PathLike, FloatPrecision
from .json import JsonLoader, JsonWriter
from core.utils.common_asserts import check_type
from core.utils.visual_iterator import VisualIterator


S = tp.TypeVar("S", bound=Structure)


class StructureFile(tp.Generic[S]):
    """
    I/O class to load and write JSON files containing efficiently stored structure data.
    Only supports ordered structures.
    """
    _data: dict[str, list[tp.Any]]
    _float_prec: FloatPrecision

    # Categorization of dict representations keys
    _dict_mandatory_keys: tuple[str, ...] = ("lattice", "sites")
    _dict_optional_keys: tuple[str, ...] = ("charge", "properties")
    _dict_keys: tuple[str, ...] = _dict_mandatory_keys + _dict_optional_keys

    # Categorization of internal storage-efficient keys
    _file_lattice_keys: tuple[str, ...] = ("a", "b", "c", "alpha", "beta", "gamma")
    _file_mandatory_site_keys: tuple[str, ...] = ("species", "site_a", "site_b", "site_c")
    _file_optional_site_keys: tuple[str, ...] = ("site_properties",)
    _file_site_keys: tuple[str, ...] = _file_mandatory_site_keys + _file_optional_site_keys
    _file_mandatory_keys: tuple[str, ...] = _file_lattice_keys + _file_mandatory_site_keys

    # NOTE: '_file_optional_keys' MUST be a dict mapping optional keys to callables returning
    # a default iterator suitable for this key for '_get_optional_key_default()' to work properly
    _file_optional_keys: dict[str, Callable[[tpe.Self], tp.Iterator[tp.Any]]] = {
        "charge": lambda self: itt.repeat(0.0, self.length), # Fastest, for immutable defaults
        "properties": lambda self: iter({} for _ in range(self.length)), # mutable defaults
        "site_properties": lambda self: iter( # site-wise attribute defaults
            [{} for _ in range(n_sites)] for n_sites in self.num_sites_per_struct_iter
        )
    }
    _file_keys: tuple[str, ...] = _file_mandatory_keys + tuple(_file_optional_keys)

    def __init__(
        self,
        data: dict[str, list[tp.Any]] | None = None,
        compressed: bool = False,
        float_prec: FloatPrecision | str | int = FloatPrecision.HIGH,
        saved_optional_keys: set[str] | None = None
    ) -> None:
        """
        I/O class to load and write JSON files containing efficiently stored structure data.
        Only supports ordered structures.

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
        # Initialize optional arguments and instance keys
        self.float_prec = float_prec

        if saved_optional_keys is None:
            saved_optional_keys = set(self._file_optional_keys)
        else:
            for key in saved_optional_keys:
                self.check_optional_key(key)

        self._saved_optional_keys = saved_optional_keys
        self._file_active_keys = set(self._file_mandatory_keys) | saved_optional_keys
        active_optional_dict_keys = {
            key for key in self._dict_optional_keys if key in saved_optional_keys
        }
        self._dict_active_keys = set(self._dict_mandatory_keys) | active_optional_dict_keys

        # No data, create new empty file
        if data is None:
            self._data = {key: [] for key in self.file_active_keys}
            return

        # Check input data conformity
        for key, val in data.items():
            check_type(key, "data key", (str,))
            check_type(val, "data value", (list,))

        # Parse and check mandatory keys only
        mandatory_keys = self._file_mandatory_keys if compressed else self._dict_mandatory_keys
        mandatory_data = _parse_dict_keys(data, mandatory_keys)
        mandatory_length = self._check_data_lengths(mandatory_data)

        if not compressed:
            # Reorganize keys for storage-efficient format
            ordered_data = (mandatory_data[k] for k in mandatory_keys)
            compressed_data_iter = (
                _compress_struct_dict(
                    {k: v for k, v in zip(mandatory_keys, struct_data, strict=True)}
                ) for struct_data in zip(*ordered_data, strict=True)
            )
            compressed_data = {key: [] for key in self._file_mandatory_keys}
            for dct in compressed_data_iter:
                for key in self._file_mandatory_keys:
                    compressed_data[key].append(dct[key])

            mandatory_data = compressed_data

        self._data = mandatory_data

        # Add saved optional keys
        for opt_key in self.saved_optional_keys:
            if opt_key in data:
                if not (opt_len:=len(data[opt_key])) == mandatory_length:
                    raise ValueError(
                        f"Saved optional key {opt_key!r} do not have expected length "
                        f"{mandatory_length}, got {opt_len}."
                    )
                self._data[opt_key] = data[opt_key]
            else:
                self._data[opt_key] = list(self._get_optional_key_default(opt_key))

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, index: int | slice) -> S | tpe.Self:
        """Get a Structure from the data at given index, or a sub-dataset from a slice."""
        if isinstance(index, slice):
            subdata = {key: self._data[key][index] for key in self.file_active_keys}
            cls = type(self)
            return cls(
                subdata,
                compressed=True,
                float_prec=self.float_prec,
                saved_optional_keys=self.saved_optional_keys
            )
        elif isinstance(index, int):
            if index < 0:
                index += self.length
            if not 0 <= index < self.length:
                raise IndexError(f"{type(self).__name__}: index out of range.")

            return self._build_structure(index)

        check_type(index, "index", (int, slice))

    def __iter__(self) -> tp.Iterator[S]:
        return iter(self._build_structure(i) for i in range(self.length))

    def __add__(self, other: tpe.Self) -> tpe.Self:
        cls = type(self)

        if not isinstance(other, cls):
            return NotImplemented

        new_dataset = {}

        # Handle data concatenation with possibly lacking optional keys
        all_active_file_keys = set(self.file_active_keys) | set(other.file_active_keys)
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

    def _build_structure(self, index: int) -> S:
        """Use data at given index to build corresponding Structure."""
        dct = _uncompress_struct_dict(
            {k: self._data[k][index] for k in self.file_active_keys}
        )
        return Structure.from_dict(dct)

    def _check_data_lengths(self, data: dict[str, list[tp.Any]] | None = None) -> int:
        """
        Check that stored data lists have same length and return it.
        If `data` is given, check lists inside `data` instead.
        """
        data = self._data if data is None else data

        if len(len_set:={len(data[key]) for key in data.keys()}) == 1:
            return len_set.pop()

        len_set_str = ", ".join(sorted(map(str, len_set)))
        raise ValueError(
            f"Some lists in data do not have the same length, got {len_set_str}."
        )

    def _check_ordering(self, structure: S) -> None:
        """Check that passed structure is ordered."""
        if not structure.is_ordered:
            raise ValueError(
                f"{type(self).__name__!r} only supports ordered structures, "
                "i.e. with only one species with occupancy 1.0 per site."
            )

    def _get_optional_key_default(self, key: str) -> tp.Iterator[tp.Any]:
        """Get an iterator containing default values for an optional data key."""
        if key in self._file_optional_keys:
            return self._file_optional_keys[key](self)
        raise ValueError(f"{key!r} is not a valid optional key.")

    def check_optional_key(self, key: str) -> None:
        """Check that the key is a valid optional key."""
        self._get_optional_key_default(key)

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
            [self._data[key] for key in self._file_lattice_keys], dtype=np.float64
        ).T
        np.round(lattice_array, decimals=self.float_prec.value, out=lattice_array)
        rounded_lists = lattice_array.T.tolist()

        for key, rounded_list in zip(self._file_lattice_keys, rounded_lists):
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

    def parse_structures(self, verbose: bool = False) -> list[S]:
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

    def add_structure(self, structure: S) -> None:
        """Add one Structure object data to the file."""
        self._check_ordering(structure)
        dct = structure.as_dict(verbosity=0)
        dct = _parse_dict_keys(dct, self.dict_active_keys)
        dct = _compress_struct_dict(dct)
        for key in self.file_active_keys:
            self._data[key].append(dct[key])

    def add_structures(self, structures: tp.Iterable[S]) -> None:
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
    def file_active_keys(self) -> set[str]:
        """Set of active internal keys containing data."""
        return self._file_active_keys

    @property
    def dict_active_keys(self) -> set[str]:
        """Set of active structure dict representation keys."""
        return self._dict_active_keys

    @property
    def length(self) -> int:
        """Number of structures contained in the file."""
        return self._check_data_lengths()

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
        params_iter = zip(*(self._data[key] for key in self._file_lattice_keys))
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
        structs_iter = zip(*(self._data[key] for key in self._file_site_keys))
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
            auto_saved_keys = set(filter(lambda key: key in data, cls._file_optional_keys))
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


def _parse_dict_keys(dct: dict[str, tp.Any], needed_keys: tp.Iterable[str]) -> dict[str, tp.Any]:
    """Verify presence of needed keys and remove other ones."""
    needed_keys = tuple(needed_keys)

    if all(keys_in_data:=[key in dct for key in needed_keys]):
        return {k: dct[k] for k in needed_keys}

    lacking_keys = ", ".join([
        needed_keys[idx]
        for idx, key_in_data in enumerate(keys_in_data) if not key_in_data
    ])
    raise KeyError(f"Some needed keys are not present: {lacking_keys}.")


def _compress_lattice_data(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Split Structure dict representation's lattice into more storage-efficient parameter keys.
    Exact reverse operation of `_uncompress_lattice_data()`.
    """
    parsed_dict = {key: dct[key] for key in dct if key != "lattice"}

    lattice = Lattice.from_dict(dct["lattice"])
    parsed_dict["a"] = lattice.a
    parsed_dict["b"] = lattice.b
    parsed_dict["c"] = lattice.c
    parsed_dict["alpha"] = lattice.alpha
    parsed_dict["beta"] = lattice.beta
    parsed_dict["gamma"] = lattice.gamma

    return parsed_dict

def _uncompress_lattice_data(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Assemble lattice parameters data into a single 'lattice' key matching
    Lattice objects dict representation. Exact reverse operation of `_compress_lattice_data()`.
    """
    parsed_dict = {key: dct[key] for key in dct if key not in StructureFile._file_lattice_keys}

    parsed_dict["lattice"] = Lattice.from_parameters(
        a=dct.pop("a"), b=dct.pop("b"), c=dct.pop("c"),
        alpha=dct.pop("alpha"), beta=dct.pop("beta"), gamma=dct.pop("gamma")
    ).as_dict()

    return parsed_dict


def _compress_sites_data(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Split internal PeriodicSite dict representations into coordinate-wise keys
    and simplify site's internal Species dict by assuming ordered site.
    Exact reverse operation of '_uncompress_sites_data()'.
    """
    parsed_dict = {key: dct[key] for key in dct if key != "sites"}

    parsed_dict["species"] = []
    parsed_dict["site_a"], parsed_dict["site_b"], parsed_dict["site_c"] = [], [], []
    parsed_dict["site_properties"] = []
    for site_dict in dct["sites"]:
        parsed_dict["species"].append(site_dict["species"][0]["element"])
        parsed_dict["site_a"].append(site_dict["abc"][0])
        parsed_dict["site_b"].append(site_dict["abc"][1])
        parsed_dict["site_c"].append(site_dict["abc"][2])
        parsed_dict["site_properties"].append(site_dict.get("properties", {}))

    return parsed_dict


def _uncompress_sites_data(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Assemble all site-wise data into a single 'sites' key matching PeriodicSite
    dict representation. Exact reverse operation of '_compress_sites_data()'.
    """
    parsed_dict = {key: dct[key] for key in dct if key not in StructureFile._file_site_keys}

    species_list = dct.pop("species")
    site_iter = zip(
        species_list,
        dct.pop("site_a"),
        dct.pop("site_b"),
        dct.pop("site_c"),
        dct.pop("site_properties", ({} for _ in range(len(species_list)))),
        strict=True
    )
    parsed_dict["sites"] = [
        {
            "species": [{"element": site_data[0], "occu": 1.0}],
            "abc": [site_data[1], site_data[2], site_data[3]],
            "properties": site_data[4]
        } for site_data in site_iter
    ]
    return parsed_dict


def _compress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into more storage-efficient organization.
    Exact reverse operation of `_uncompress_struct_dict()`.
    """
    file_dict = _compress_lattice_data(dct)
    file_dict = _compress_sites_data(file_dict)
    return file_dict


def _uncompress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into object reconstruction-friendly organization.
    Exact reverse operation of `_compress_struct_dict()`.
    """
    struct_dict = _uncompress_lattice_data(dct)
    struct_dict = _uncompress_sites_data(struct_dict)
    return struct_dict
