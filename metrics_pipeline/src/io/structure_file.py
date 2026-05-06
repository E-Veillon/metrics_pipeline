"""
Class handling specific JSON files containing structures data in a predictable and compact
way to avoid ambiguities in data communications between pipeline steps.
"""

import typing as tp
import typing_extensions as tpe
from collections.abc import Mapping
import functools as ft
import itertools as itt

from pymatgen.core import Structure, Lattice, PeriodicSite

from .io_base import PathLike
from .json import JsonLoader, JsonWriter
from src.utils.common_asserts import check_type


class StructureFile:
    """
    I/O class to load and write JSON files containing efficiently stored structure data.
    """
    _struct_file_keys = (
        "a", "b", "c", "alpha", "beta", "gamma", "charge", "properties",
        "species", "site_a", "site_b", "site_c", "site_properties"
    )
    _struct_dict_keys = (
        "lattice", "charge", "properties", "sites"
    )
    _data: dict[str, list[tp.Any]]
    _float_prec: tp.Literal["low", "medium", "high"]

    def __init__(
        self,
        data: dict[str, list[tp.Any]] | None = None,
        compressed: bool = False,
        float_prec: tp.Literal["low", "medium", "high"] = "high",
        save_charge: bool = True,
        save_properties: bool = True,
        save_site_properties: bool = True,
        _from_file: bool = False
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

        float_prec: 'low', 'medium', 'high'
            What level of precision to apply to lattice lengths and angles and site positions
            values when storing structure data into the file:
            - 'low' rounds to 4 decimal places
            (0.01 pm for lengths and positions, 1e-4 degree for angles);
            - 'medium' rounds to 8 decimal places;
            - 'high' rounds to 16 decimal places (numerical precision level).

            WARNING: Lost precision is irreversible and will NOT be recovered when retrieving data.
            Defaults to 'high' to avoid unintended precision loss.

        save_charge: bool
            Whether to save the overall charge of structures in the file. Set to False if storage
            efficiency is critical (e.g. lot of structures) and this optional data is not useful
            to you. Defaults to True.

        save_properties: bool
            Whether to save the structures properties dicts in the file. Set to False if storage
            efficiency is critical (e.g. lot of structures) and this optional data is not useful
            to you. Defaults to True.

        save_site_properties: bool
            Whether to save the structures site properties dicts in the file. Set to False if
            storage efficiency is critical (e.g. lot of structures) and this optional data is
            not useful to you. Defaults to True.

        _from_file: bool
            Internal flag to indicate that passed data was loaded from a file. Do not set manually.
        """
        self.float_prec = float_prec
        self._save_charge = save_charge
        self._save_properties = save_properties
        self._save_site_properties = save_site_properties

        if data is None:
            self._data = {}
            for key in self.active_file_keys:
                self._data[key] = []
            return

        keys_to_check = self.active_file_keys if compressed else self._struct_dict_keys
        _check_dict_keys(data, keys_to_check)

        if not len(len_set:={len(data[key]) for key in keys_to_check}) == 1:
            len_set_str = ", ".join(sorted(map(str, len_set)))
            raise ValueError(
                f"Some iterables in data do not have the same length, got {len_set_str}."
            )
        data = data if compressed else self._reorganize_data(data, algo="compress")

        if not save_charge:
            data.pop("charge", None)
        
        if not save_properties:
            data.pop("properties", None)
        
        if not save_site_properties:
            data.pop("site_properties", None)

        self._data = data

        if _from_file:
            # If a file was loaded without some optional keys but we now want to add
            # new data with these optional data saved, generate default values for
            # previously unsaved optional data to maintain file consistency
            if save_charge and "charge" not in data:
                self._data["charge"] = list(
                    self._get_optional_key_default("charge")
                )
            if save_properties and "properties" not in data:
                self._data["properties"] = list(
                    self._get_optional_key_default("properties")
                )
            if save_site_properties and "site_properties" not in data:
                self._data["site_properties"] = list(
                    self._get_optional_key_default("site_properties")
                )

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, index: int | slice) -> Structure | tpe.Self:
        """Get a Structure from the data at given index."""
        if isinstance(index, slice):
            subdata = {key: self._data[key][index] for key in self.active_file_keys}
            cls = type(self)
            return cls(
                subdata,
                compressed=True,
                save_charge=self._save_charge,
                save_properties=self._save_properties,
                save_site_properties=self._save_site_properties
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
        if not isinstance(other, type(self)):
            return NotImplemented

        # Handle optional keys
        diff_keys = (
            self._save_charge - other._save_charge,
            self._save_properties - other._save_properties,
            self._save_site_properties - other._save_site_properties
        )

        for idx, key in enumerate(("charge", "properties", "site_properties")):
            if diff_keys[idx] == 0: # Both have or don't have it
                continue
            if diff_keys[idx] == -1: # self does not have it but other does
                self._data[key] = list(self._get_optional_key_default(key))
                setattr(self, f"_save_{key}", True)
            else: # self has it but not other
                other._data[key] = list(other._get_optional_key_default(key))
                setattr(other, f"_save_{key}", True)

        # Safety check of keys (should never trigger)
        if self.active_file_keys != other.active_file_keys:
            self_supp_keys = ', '.join(sorted(
                set(self.active_file_keys) - set(other.active_file_keys)
            ))
            other_supp_keys = ', '.join(sorted(
                set(other.active_file_keys) - set(self.active_file_keys)
            ))
            raise AssertionError(
                "Somehow databases have incompatible keys after processing optional keys, "
                "this should not happen and is likely a bug. See problematic keys below:\n"
                f"- left operand: {self_supp_keys}\n"
                f"- right operand: {other_supp_keys}\n"
            )

        # Concatenate each data list
        new_data = {
            key: self._data[key] + other._data[key]
            for key in self.active_file_keys
        }
        cls = type(self)
        return cls(
            new_data,
            compressed=True,
            save_charge=self._save_charge,
            save_properties=self._save_properties,
            save_site_properties=self._save_site_properties
        )

    def _reorganize_data(
        self, data: Mapping[str, tp.Iterable[tp.Any]], *, algo: tp.Literal["compress", "uncompress"]
    ) -> dict[str, list[tp.Any]]:
        """
        Reorganize data between storage-efficient and object reconstruction-friendly formats.
        """
        if algo == "compress":
            parser_fn = ft.partial(_compress_struct_dict, float_prec=self.float_prec)
            init_key_list = self._struct_dict_keys
            target_key_list = self.active_file_keys
        elif algo == "uncompress":
            parser_fn = _uncompress_struct_dict
            init_key_list = self.active_file_keys
            target_key_list = self._struct_dict_keys
        else:
            raise ValueError(f"'algo' must be either 'compress' or 'uncompress', got {algo!r}.")

        ordered_data = (data[k] for k in init_key_list)

        try:
            parsed_data = [
                parser_fn({k: v for k, v in zip(init_key_list, struct_data, strict=True)})
                for struct_data in zip(*ordered_data, strict=True)
            ]
        except ValueError as exc:
            raise ValueError(
                "Data reorganization failed: some iterables in the data dict "
                "do not have the same length."
            ) from exc

        reorganized_data: dict[str, list[tp.Any]] = {}
        for key in target_key_list:
            reorganized_data[key] = [struct_dict[key] for struct_dict in parsed_data]

        return reorganized_data

    def _build_structure(self, index: int) -> Structure:
        """Use data at given index to build corresponding Structure."""
        dct = self._reorganize_data(
            {k: self._data[k][index:index+1] for k in self.active_file_keys},
            algo="uncompress"
        )
        dct = {k: v[0] for k, v in dct.items()}
        return Structure.from_dict(dct)

    def _check_ordering(self, structure: Structure) -> None:
        """Check that passed structure is ordered."""
        if not structure.is_ordered:
            raise ValueError(
                f"{type(self).__name__!r} only supports ordered structures, "
                "i.e. with only one species with occupancy 1.0 per site."
            )

    def _get_optional_key_default(self, key: str) -> tp.Iterator[tp.Any]:
        match key:
            case "charge":
                return itt.repeat(0.0, self.length)
            case "properties":
                return iter({} for _ in range(self.length))
            case "site_properties":
                return iter(
                    [{} for _ in range(n_sites)]
                    for n_sites in self.num_sites_per_struct_iter
                )
            case x:
                raise ValueError(f"{x!r} is not a valid optional key.")

    def parse_structures(self) -> list[Structure]:
        """Build all `Structure` objects from file data, in order."""
        return list(self)

    def add_structure(self, structure: Structure) -> None:
        """Add one Structure object data to the file."""
        self._check_ordering(structure)
        dct = structure.as_dict(verbosity=0)
        dct = _compress_struct_dict(dct, float_prec=self.float_prec)
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
    def active_file_keys(self) -> tuple[str, ...]:
        """Get internal keys where data is actively saved."""
        file_keys = list(self._struct_file_keys)
        if not self._save_charge:
            file_keys.remove("charge")
        if not self._save_properties:
            file_keys.remove("properties")
        if not self._save_site_properties:
            file_keys.remove("site_properties")
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
        params_iter = zip(*(self._data[key] for key in ("a", "b", "c", "alpha", "beta", "gamma")))
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
        if self._save_charge:
            return iter(self._data["charge"])
        return self._get_optional_key_default("charge")

    @property
    def all_properties_iter(self) -> tp.Iterator[dict]:
        """
        Lazy iterator of all stored structure properties dicts in order. If `save_properties`
        was set to False, builds an iterator of empty properties dicts coherent with file data
        for file iteration consistency.
        """
        if self._save_properties:
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
        if self._save_site_properties:
            return iter(self._data["site_properties"])
        return self._get_optional_key_default("site_properties")

    @property
    def all_sites_iter(self) -> tp.Iterator[list[PeriodicSite]]:
        """
        Lazy iterator of all lists of sites (one list per structure)
        rebuilt from data, in order.
        """
        sites_keys = ("species", "site_a", "site_b", "site_c", "site_properties")
        structs_iter = zip(*(self._data[key] for key in sites_keys))
        lattices_iter = self.all_lattices_iter

        for sites_data_tuple, lattice in zip(structs_iter, lattices_iter):
            sites_iter = zip(*sites_data_tuple)
            yield [
                PeriodicSite(
                    elt, [site_a, site_b, site_c], lattice, site_props
                ) for elt, site_a, site_b, site_c, site_props in sites_iter
            ]

    @property
    def float_prec(self) -> tp.Literal["low", "medium", "high"]:
        """Initialized compressed float precision level."""
        return self._float_prec

    @float_prec.setter
    def float_prec(self, value: tp.Literal["low", "medium", "high"]) -> None:
        """Set a new float precision level. Only applies to further added data."""
        match value:
            case "low" | "medium" | "high":
                self._float_prec = value
            case _:
                raise ValueError(
                    "'float_prec' only supports 'low', 'medium' and 'high', "
                    f"got {value!r}."
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
            Whether to pass optional data args automatically to match keys that are
            present in file. Defaults to True. Avoids generating default empty data
            at initialization by default whenever loading a file. Overrides any passed
            optional data argument. Set to False to pass your own combination of optional
            data as keyword arguments.

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

        if match_optional_keys:
            kwargs["save_charge"] = "charge" in data
            kwargs["save_properties"] = "properties" in data
            kwargs["save_site_properties"] = "site_properties" in data

        return cls(data, compressed=True, _from_file=True, **kwargs)

    def write_file(self, filepath: PathLike) -> None:
        """Write stored data to a compact JSON file."""
        written_data = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "data": self._data
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


def _compress_struct_dict(
    dct: dict[str, tp.Any], float_prec: tp.Literal["low", "medium", "high"] = "high"
) -> dict[str, tp.Any]:
    """
    Reorganize structure data into more storage-efficient organization.
    Exact reverse operation of `_uncompress_struct_dict()`.

    Parameters
    ----------
    dct: dict[str, Any]
        Object reconstruction-friendly dict to compress to storage-efficient format.

    float_prec: 'low', 'medium', 'high'
        What level of precision to apply to lattice lengths and angles and site positions values.
        - 'low' rounds to 4 decimal places (0.01 pm for lengths and positions, 1e-4 degree for angles);
        - 'medium' rounds to 8 decimal places;
        - 'high' rounds to 16 decimal places (numerical precision level).

        WARNING: Lost precision is irreversible and will NOT be recovered by decompression.
        Defaults to 'high' to avoid unintended precision loss.

    Returns
    -------
    dcit[str, tp.Any]
        Storage-efficient formatted structure dict.
    """
    match float_prec:
        case "low": ndigits = 4
        case "medium": ndigits = 8
        case "high": ndigits = 16
        case _:
            raise ValueError(
                f"'float_prec' only supports 'low', 'medium' and 'high', got {float_prec!r}."
            )
    rounder = ft.partial(round, ndigits=ndigits)

    file_dict = {}
    file_dict["charge"] = dct.pop("charge")
    file_dict["properties"] = dct.pop("properties")

    # Flatten lattice parameters
    lattice = Lattice.from_dict(dct["lattice"])
    file_dict["a"] = rounder(lattice.a)
    file_dict["b"] = rounder(lattice.b)
    file_dict["c"] = rounder(lattice.c)
    file_dict["alpha"] = rounder(lattice.alpha)
    file_dict["beta"] = rounder(lattice.beta)
    file_dict["gamma"] = rounder(lattice.gamma)

    # Flatten sites data
    file_dict["species"] = []
    file_dict["site_a"], file_dict["site_b"], file_dict["site_c"] = [], [], []
    file_dict["site_properties"] = []
    for site_dict in dct["sites"]:
        file_dict["species"].append(site_dict["species"][0]["element"])
        file_dict["site_a"].append(rounder(site_dict["abc"][0]))
        file_dict["site_b"].append(rounder(site_dict["abc"][1]))
        file_dict["site_c"].append(rounder(site_dict["abc"][2]))
        file_dict["site_properties"].append(site_dict["properties"])

    return file_dict


def _uncompress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into object reconstruction-friendly organization.
    Exact reverse operation of `_compress_struct_dict()` (except for float precision compression).
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
