"""
Class handling specific JSON files containing GenMat structures data in a predictable and compact
way to avoid ambiguities in data communications between pipeline steps.
"""

import typing as tp
import typing_extensions as tpe
from collections.abc import Mapping
import functools as ft

from pymatgen.core import Structure, Lattice

from .io_base import PathLike
from .json import JsonLoader, JsonWriter
from src.utils.common_asserts import check_type
from src.utils.genmat_data import GenMatStructure
from src.utils.spg_data import Spacegroup


class GenMatFile:
    """
    I/O class to load and write JSON files containing structure data for the GenMat library.
    """
    _genmat_file_keys = (
        "name", "special_keys", "spacegroup", "energy", "energy_above_hull",
        "a", "b", "c", "alpha", "beta", "gamma", "charge", "properties",
        "species", "site_a", "site_b", "site_c", "site_properties"
    )
    _genmat_struct_keys = (
        "name", "special_keys", "spacegroup", "energy", "energy_above_hull", "structure"
    )
    _data: dict[str, list[tp.Any]]

    def __init__(
        self,
        data: dict[str, list[tp.Any]] | None = None,
        compressed: bool = False,
        float_prec: tp.Literal["low", "medium", "high"] = "high"
    ) -> None:
        """
        I/O class to load and write JSON files containing structure data for the GenMat library.

        Parameters
        ----------
        data: dict[str, list[Any]], optional
            GenMatStructrure data to put into the file. Must be a dict with the same keys
            as a dict generated with `GenMatStrtucture.as_dict()`, containing lists of equal
            lengths where the same index contains data for the same structure. If not given,
            instantiate an empty file that can be filled in via `add_structure()` and
            `add_structures()` methods.

        compressed: bool
            Whether given data already follows the storage-efficient dict formatting, which
            optimizes the `structure` attribute saved data into several easier-to-store values.
            Defaults to False.

        float_prec: 'low', 'medium', 'high'
            What level of precision to apply to lattice lengths and angles and site positions values
            when storing structure data into the file:
            - 'low' rounds to 4 decimal places (0.01 pm for lengths and positions, 1e-4 degree for angles);
            - 'medium' rounds to 8 decimal places;
            - 'high' rounds to 16 decimal places (numerical precision level).

            WARNING: Lost precision is irreversible and will NOT be recovered when retrieving data.
            Defaults to 'high' to avoid unintended precision loss.
        """
        self.float_prec = float_prec
        if data is None:
            self._data = {}
            for key in self._genmat_file_keys:
                self._data[key] = []
            return

        keys_to_check = self._genmat_file_keys if compressed else self._genmat_struct_keys
        _check_dict_keys(data, keys_to_check)

        if not len(len_set:={len(data[key]) for key in keys_to_check}) == 1:
            len_set_str = ", ".join(sorted(map(str, len_set)))
            raise ValueError(
                f"Some iterables in data do not have the same length, got {len_set_str}."
            )

        self._data = data if compressed else self._reorganize_data(data, algo="compress")

    def __len__(self) -> int:
        return len(self.all_names)

    def __getitem__(self, index: int) -> GenMatStructure:
        """Get a GenMatStructure from the data at given index."""
        check_type(index, "index", (int,))
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(f"{type(self).__name__}: index out of range.")

        return self._build_structure(index)

    def __iter__(self) -> tp.Iterator[GenMatStructure]:
        return iter(self._build_structure(i) for i in range(len(self)))

    def _reorganize_data(
        self, data: Mapping[str, tp.Iterable[tp.Any]], *, algo: tp.Literal["compress", "uncompress"]
    ) -> dict[str, list[tp.Any]]:
        """
        Reorganize data between storage-efficient and object reconstruction-friendly formats.
        """
        if algo == "compress":
            parser_fn = partial(_compress_struct_dict, float_prec=self.float_prec)
            init_key_list = self._genmat_struct_keys
            target_key_list = self._genmat_file_keys
        elif algo == "uncompress":
            parser_fn = _uncompress_struct_dict
            init_key_list = self._genmat_file_keys
            target_key_list = self._genmat_struct_keys
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

    def _build_structure(self, index: int) -> GenMatStructure:
        """Use data at given index to build corresponding GenMatStructure."""
        dct = self._reorganize_data(
            {k: self._data[k][index:index+1] for k in self._genmat_file_keys},
            algo="uncompress"
        )
        dct = {k: v[0] for k, v in dct.items()}
        return GenMatStructure.from_dict(dct)

    def _check_ordering(self, structure: GenMatStructure) -> None:
        """Check that passed structure is ordered."""
        if not structure.is_ordered:
            raise ValueError(
                f"{type(self).__name__!r} only supports ordered structures, "
                "i.e. with only one species with occupancy 1.0 per site. "
                f"Disordered structure {structure.name!r} is not supported."
            )

    def parse_structures(self) -> list[GenMatStructure]:
        """Build all `GenMatStructure` objects from file data, in order."""
        return list(self)

    def parse_pmg_structures(self) -> list[Structure]:
        """Build all standard pymatgen `Structure` objects from file data, in order."""
        return [gstruct.as_structure() for gstruct in self]

    def add_structure(self, structure: GenMatStructure) -> None:
        """Add a `GenMatStructure` data to the file."""
        self._check_ordering(structure)
        dct = structure.as_dict(verbosity=0)
        dct = _compress_struct_dict(dct, float_prec=self.float_prec)
        _check_dict_keys(dct, self._genmat_file_keys)
        for key in self._genmat_file_keys:
            self._data[key].append(dct[key])

    def add_structures(self, structures: tp.Iterable[GenMatStructure]) -> None:
        """
        Convenient method to add a full iterable of `GenMatStructure` data to the file.
        Order is preserved.
        """
        for structure in structures:
            self.add_structure(structure)

    @property
    def float_prec(self) -> str:
        """Initialized compressed float precision level."""
        return self._float_prec

    @float_prec.setter
    def float_prec(self, value) -> None:
        """Set a new float precision level. Only applies to further added data."""
        match value:
            case "low" | "medium" | "high":
                self._float_prec = value
            case _:
                raise ValueError(
                    f"'float_prec' only supports 'low', 'medium' and 'high', got {float_prec!r}."
                )

    @property
    def all_names(self) -> list[str]:
        """List of all stored structure names."""
        return self._data["name"]

    @property
    def all_special_keys(self) -> list[dict[str, tp.Any]]:
        """List of all stored special keys dictionaries."""
        return self._data["special_keys"]

    @property
    def all_spacegroups(self) -> list[Spacegroup]:
        """List of all stored spacegroups, converted into Spacegroup objects on the fly."""
        return [Spacegroup.from_dict(dct) for dct in self._data["spacegroup"]]

    @property
    def all_energies(self) -> list[float | None]:
        """List of all stored energy values."""
        return self._data["energy"]

    @property
    def all_energies_above_hull(self) -> list[float | None]:
        """List of all stored energy above hull values."""
        return self._data["energy_above_hull"]

    @classmethod
    def from_file(cls, filepath: PathLike) -> tpe.Self:
        """Load a JSON file of GenMat data."""
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
        return cls(data, compressed=True)

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
    Reorganize structure data into more storage-efficient organization in place.
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

    struct_dict = dct.pop("structure")
    dct["charge"] = struct_dict["charge"]
    dct["properties"] = struct_dict["properties"]

    # Flatten lattice parameters
    lattice = Lattice.from_dict(struct_dict["lattice"])
    dct["a"] = rounder(lattice.a)
    dct["b"] = rounder(lattice.b)
    dct["c"] = rounder(lattice.c)
    dct["alpha"] = rounder(lattice.alpha)
    dct["beta"] = rounder(lattice.beta)
    dct["gamma"] = rounder(lattice.gamma)

    # Flatten sites data
    dct["species"] = []
    dct["site_a"], dct["site_b"], dct["site_c"] = [], [], []
    dct["site_properties"] = []
    for site_dict in struct_dict["sites"]:
        dct["species"].append(site_dict["species"][0]["element"])
        dct["site_a"].append(rounder(site_dict["abc"][0]))
        dct["site_b"].append(rounder(site_dict["abc"][1]))
        dct["site_c"].append(rounder(site_dict["abc"][2]))
        dct["site_properties"].append(site_dict["properties"])

    return dct


def _uncompress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into object reconstruction-friendly organization in place.
    Exact reverse operation of `_compress_struct_dict()` (except for float precision compression).
    """
    struct_dict = {}
    struct_dict["charge"] = dct.pop("charge")
    struct_dict["properties"] = dct.pop("properties")

    # Rebuild pymatgen-style lattice dict
    struct_dict["lattice"] = Lattice.from_parameters(
        a=dct.pop("a"), b=dct.pop("b"), c=dct.pop("c"),
        alpha=dct.pop("alpha"), beta=dct.pop("beta"), gamma=dct.pop("gamma")
    ).as_dict()

    # Rebuild pymatgen-style site dict
    site_iter = zip(
        dct.pop("species"),
        dct.pop("site_a"),
        dct.pop("site_b"),
        dct.pop("site_c"),
        dct.pop("site_properties"),
        strict=True
    )
    struct_dict["sites"] = [
        {
            "species": [{"element": site_data[0], "occu": 1.0}],
            "abc": [site_data[1], site_data[2], site_data[3]],
            "properties": site_data[4]
        } for site_data in site_iter
    ]

    dct["structure"] = struct_dict

    return dct
