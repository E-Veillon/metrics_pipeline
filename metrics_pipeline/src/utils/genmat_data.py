"""GenMat pipeline specific functionnalities."""

import re
import typing as tp
import typing_extensions as tpe
from collections.abc import Mapping

from pymatgen.core.structure import Structure, Composition, Lattice
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .spg_data import Spacegroup
from .visual_iterator import VisualIterator
from .common_asserts import check_type
from src.computations.local import get_density
from src.io import JsonLoader, JsonWriter, PathLike


_ELT_WITH_INDEX_REGEX = r"[A-Z][a-z]?\d*"
_GENMAT_NAME_PATTERN = re.compile(
    fr"^(\d+)_(?:{_ELT_WITH_INDEX_REGEX}|\((?:{_ELT_WITH_INDEX_REGEX})+\)\d*)+$"
)


def is_genmat_name(name: str) -> bool:
    """
    Whether given name follows GenMat naming convention.
    """
    return isinstance(name, str) and _GENMAT_NAME_PATTERN.fullmatch(name) is not None


def check_genmat_name(name: str) -> None:
    """
    Verify that given name follows GenMat naming convention.
    Raises a `ValueError` if it does not.
    """
    if is_genmat_name(name):
        return
    raise ValueError(
        f"{name!r} does not follow GenMat naming conventions. "
        "Please make sure your data was preprocessed with the 'preprocess.py' "
        "script before going further."
    )


def get_genmat_name_idx(name: str) -> int:
    """
    Extract the index part of a GenMat name.
    """
    check_genmat_name(name)
    match = _GENMAT_NAME_PATTERN.match(name)
    if match is None:
        raise ValueError(f"Failed to extract index from {name!r}")
    return int(match.group(1))


def generate_genmat_names(structures: list[Structure]) -> list[Structure]:
    """
    Generate GenMat indexed names based on the position of each structure in the list.
    The new name is stored inside the 'header' key of structures 'properties' attribute.
    If another name already exists in 'header', it is moved to the '_original_header' key.
    """
    iterator: VisualIterator[tuple[int, Structure]] = VisualIterator.from_big_iterator(
        enumerate(structures), n_elts=len(structures),
        desc="Generating GenMat names", unit="generated", percent=True
    )
    for idx, structure in iterator:
        if old_header:=structure.properties.get('header'):
            structure.properties["_original_header"] = old_header
        structure.properties["header"] = f"{idx}_{structure.reduced_formula}"

    return structures


class GenMatPDEntry(PDEntry):
    """
    Simple extension of pymatgen's `PDEntry` objects making their `name` attribute mandatory
    and must follow GenMat naming convention, and add an `energy_above_hull` property to
    store computed energy above hull in eV/atom when compared with a `PhaseDiagram` object.

    Attributes
    ----------
    composition: Composition
        The composition associated with the PDEntry.
    energy: float
        The energy associated with the entry.
    name: str
        A name for the entry. This is the string shown in the phase diagrams.
        Must follow GenMat naming convention.
    energy_above_hull: float | None
        Computed energy above hull of the entry when compared to a `PhaseDiagram` object.
    attribute: MSONable
        A arbitrary attribute. Can be used to specify that the
        entry is a newly found compound, or to specify a particular label for
        the entry, etc. An attribute can be anything but must be MSONable.
    """
    def __init__(
        self,
        composition: Composition,
        energy: float,
        name: str,
        energy_above_hull: float | None = None,
        attribute: object = None
    ):
        """
        Simple extension of pymatgen's `PDEntry` objects making their `name` attribute
        mandatory and must follow GenMat naming convention, and add an `energy_above_hull`
        property to store computed energy above hull in eV/atom when compared with a
        `PhaseDiagram` object.

        Parameters
        ----------
        composition: Composition
            The composition associated with the PDEntry.

        energy: float
            The energy associated with the entry.

        name: str
            A name for the entry. This is the string shown in the phase diagrams.
            Must follow GenMat naming convention.

        energy_above_hull: float, optional
            Computed energy above hull of the entry when compared to a `PhaseDiagram` object.

        attribute: MSONable, optional
            A arbitrary attribute. Can be used to specify that the
            entry is a newly found compound, or to specify a particular label for
            the entry, etc. An attribute can be anything but must be MSONable.
        """
        check_genmat_name(name)
        super().__init__(composition, energy, name, attribute)
        check_type(energy_above_hull, "energy_above_hull", (float, type(None)))
        self._energy_above_hull = energy_above_hull

    @property
    def energy_above_hull(self) -> float | None:
        return self._energy_above_hull

    @energy_above_hull.setter
    def energy_above_hull(self, value: float | None) -> None:
        check_type(value, "energy_above_hull", (float, type(None)))
        self._energy_above_hull = value

    @energy_above_hull.deleter
    def energy_above_hull(self) -> None:
        self._energy_above_hull = None

    def compute_energy_above_hull(self, pdiagram: PhaseDiagram, **kwargs: tp.Any) -> None:
        """
        Compute energy above hull of the entry relative to a phase diagram in eV/atom.
        The result is saved in the `energy_above_hull` attribute.

        Parameters
        ----------
        pdiagram: PhaseDiagram
            The phase diagram to compare the entry to.

        kwargs: Any
            Additional keyword arguments to pass to `PhaseDiagram.get_decomp_and_e_above_hull()`.
        """
        kwargs.setdefault("allow_negative", True)
        kwargs.setdefault("check_stable", False)
        e_above_hull = pdiagram.get_e_above_hull(self, **kwargs)
        self.energy_above_hull = e_above_hull


class GenMatStructure(Structure):
    """
    An extension of pymatgen's `Structure` objects to store GenMat name,
    optional saved CIF labels, and useful properties for metrics computations.
    Only extension attributes are shown below, but the class also inherits all
    `Structure` attributes.

    Attributes
    ----------
    name: str
        The name of the structure, following GenMat naming convention.

    special_keys: dict[str, Any]
        A dictionary to store any saved CIF labels.

    spacegroup: Spacegroup
        The spacegroup of the structure.

    energy: float | None
        The energy of the structure, in eV.

    energy_per_atom: float | None
        The energy per atom of the structure, in eV/atom.

    energy_above_hull: float | None
        The energy above hull of the structure, in eV/atom.

    density: float
        The volumic mass of the structure, in g/cm³.
    """
    def __init__(
        self,
        name: str,
        special_keys: dict[str, tp.Any] | None = None,
        spacegroup: Spacegroup | None = None,
        energy: float | None = None,
        energy_above_hull: float | None = None,
        **kwargs
    ) -> None:
        """
        An extension of pymatgen's `Structure` objects to store GenMat name,
        optional saved CIF labels, and useful properties for metrics computations.

        Parameters
        ----------
        name: str
            The name of the structure, following GenMat naming convention.

        special_keys: dict[str, Any], optional
            A dictionary to store any saved CIF labels.

        spacegroup: Spacegroup, optional
            The spacegroup of the structure.

        energy: float, optional
            The energy of the structure, in eV.

        energy_above_hull: float, optional
            The energy above hull of the structure, in eV/atom.

        kwargs: Any
            Keyword arguments to pass to the `Structure` constructor.

        Notes
        -----
        It is likely one does not want to build structures from scratch, but rather
        convert existing `Structure` objects to `GenMatStructure` for metrics computations.
        See the `from_structure()` method for seemless conversions safe of any structure
        data modification.
        """
        check_genmat_name(name)
        check_type(special_keys, "special_keys", (dict, type(None)))
        check_type(spacegroup, "spacegroup", (Spacegroup, type(None)))
        check_type(energy, "energy", (float, type(None)))
        check_type(energy_above_hull, "energy_above_hull", (float, type(None)))
        self.name = name
        self.special_keys = special_keys or {}
        self.spacegroup = spacegroup or Spacegroup(0)
        self.energy = energy
        self.energy_above_hull = energy_above_hull

        super().__init__(**kwargs)

    @property
    def name_index(self) -> int:
        """The index part of the structure name."""
        match = _GENMAT_NAME_PATTERN.match(self.name)
        if match is None:
            raise RuntimeError(
                f"Structure name {self.name!r} failed index matching. "
                "This should never happen as names are validated at initialization."
            )
        return int(match.group(1))

    @property
    def entry(self) -> GenMatPDEntry:
        """The `GenMatPDEntry` object associated with the structure."""
        if self.energy is None:
            raise ValueError(
                f"Energy must be defined to generate entry for structure {self.name!r}."
            )
        return GenMatPDEntry(self.composition, self.energy, self.name, self.energy_above_hull)

    @property
    def energy_per_atom(self) -> float | None:
        """The energy per atom of the structure, in eV/atom."""
        if self.energy is None:
            return None
        return self.energy / self.composition.num_atoms

    @property
    def density(self) -> float:
        """Volumic mass of the structure in g.cm⁻³."""
        return get_density(self)

    def compute_energy_above_hull(self, pdiagram: PhaseDiagram, **kwargs: tp.Any) -> None:
        """
        Compute energy above hull of the structure relative to a phase diagram in eV/atom.
        The result is saved in the `energy_above_hull` attribute.

        Parameters
        ----------
        pdiagram: PhaseDiagram
            The phase diagram to compare the structure to.

        kwargs: Any
            Additional keyword arguments to pass to `PhaseDiagram.get_decomp_and_e_above_hull()`.
        """
        computed_entry = self.entry
        computed_entry.compute_energy_above_hull(pdiagram, **kwargs)
        self.energy_above_hull = computed_entry.energy_above_hull

    def as_dict(self, **kwargs) -> dict[str, tp.Any]:
        """
        Get a dictionary representation of the GenMatStructure object.
        Additional keyword arguments are directly passed to `Structure.as_dict()`.
        """
        struct_dct = super().as_dict(**kwargs)
        struct_dct.pop("@module", None)
        struct_dct.pop("@class", None)
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "name": self.name,
            "special_keys": self.special_keys,
            "spacegroup": self.spacegroup.as_dict(),
            "energy": self.energy,
            "energy_above_hull": self.energy_above_hull,
            "structure": struct_dct
        }
        return dct

    def as_structure(self) -> Structure:
        """
        Rebuild a standard pymatgen `Structure` from the GenMatStructure.
        All data related to the `Structure` object are preserved.
        """
        return Structure(
            lattice=self.lattice,
            species=self.species,
            coords=self.frac_coords,
            charge=self.charge,
            site_properties=self.site_properties,
            labels=self.labels,
            properties=self.properties
        )

    @classmethod
    def from_dict(cls, dct: dict[str, tp.Any]) -> tpe.Self:
        """Create a GenMatStructure object from a dictionary representation."""
        return cls(
            name=dct["name"],
            special_keys=dct.get("special_keys", {}),
            spacegroup=Spacegroup.from_dict(dct["spacegroup"]),
            energy=dct.get("energy"),
            energy_above_hull=dct.get("energy_above_hull"),
            **dct["structure"]
        )

    @classmethod
    def from_structure(cls, name: str, structure: Structure, **kwargs) -> tpe.Self:
        """
        Build a `GenMatStructure` from an existing `Structure` object.
        This conversion is only meant as an extension of the original structure to
        the features of `GenMatStructure` class. The `Structure` data will be kept as-is
        without any additional check or transformation purposefully for data fidelity.

        Parameters
        ----------
        name: str
            The name of the structure, following GenMat naming convention.

        structure: Structure
            `Structure` object to convert into a `GenMatStructure`.

        kwargs: Any
            Additional keyword arguments to pass to the constructor.
        """
        special_keys = kwargs.pop("special_keys", None)
        spacegroup = kwargs.pop("spacegroup", None)
        energy = kwargs.pop("energy", None)
        energy_above_hull = kwargs.pop("energy_above_hull", None)
        return cls(
            name, special_keys, spacegroup, energy, energy_above_hull,
            lattice = structure.lattice,
            species = structure.species,
            coords = structure.frac_coords,
            charge = structure.charge,
            site_properties = structure.site_properties,
            labels = structure.labels,
            properties = structure.properties
        )


PDEntryLike: tp.TypeAlias = PDEntry | GenMatPDEntry
StructureLike: tp.TypeAlias = Structure | GenMatStructure


def generate_genmat_structures(
    structures: list[Structure], special_keys: list[str] | None = None
) -> list[GenMatStructure]:
    """
    Generate GenMatStructure objects from a list of Structure objects.
    The GenMat name is generated based on the position of each structure in the list,
    and stored in the 'name' attribute of the GenMatStructure object.
    """
    if special_keys is None:
        special_keys = []

    genmat_structures = [
        GenMatStructure.from_structure(
            name=f"{idx}_{structure.reduced_formula}",
            structure=structure,
            special_keys={k: structure.properties.get(k) for k in special_keys}
        ) for idx, structure in enumerate(structures)
    ]
    return genmat_structures


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

    @classmethod
    def _reorganize_data(
        cls, data: Mapping[str, tp.Iterable[tp.Any]], *, algo: tp.Literal["compress", "uncompress"]
    ) -> dict[str, list[tp.Any]]:
        """
        Reorganize data between storage-efficient and object reconstruction-friendly formats.
        """
        if algo == "compress":
            parser_fn = _compress_struct_dict
            init_key_list = cls._genmat_struct_keys
            target_key_list = cls._genmat_file_keys
        elif algo == "uncompress":
            parser_fn = _uncompress_struct_dict
            init_key_list = cls._genmat_file_keys
            target_key_list = cls._genmat_struct_keys
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

    def __init__(self, data: dict[str, list[tp.Any]] | None = None, compressed: bool = False) -> None:
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
        """
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
        """Build all `GenMatStructure` objects from file data."""
        return list(self)

    def parse_pmg_structures(self) -> list[Structure]:
        """Build all standard pymatgen `Structure` objects from file data."""
        return [gstruct.as_structure() for gstruct in self]

    def add_structure(self, structure: GenMatStructure) -> None:
        """Add a `GenMatStructure` data to the file."""
        self._check_ordering(structure)
        dct = structure.as_dict(verbosity=0)
        compressed_dct = _compress_struct_dict(dct)
        _check_dict_keys(compressed_dct, self._genmat_file_keys)
        for key in self._genmat_file_keys:
            self._data[key].append(compressed_dct[key])

    def add_structures(self, structures: tp.Iterable[GenMatStructure]) -> None:
        """
        Convenient method to add a full iterable of `GenMatStructure` data to the file.
        Order is preserved.
        """
        for structure in structures:
            self.add_structure(structure)

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


def _compress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into more storage-efficient organization in place.
    Exact reverse operation of `_uncompress_struct_dict()`.
    """
    struct_dict = dct.pop("structure")
    dct["charge"] = struct_dict["charge"]
    dct["properties"] = struct_dict["properties"]

    # Flatten lattice parameters
    lattice = Lattice.from_dict(struct_dict["lattice"])
    dct["a"] = lattice.a
    dct["b"] = lattice.b
    dct["c"] = lattice.c
    dct["alpha"] = lattice.alpha
    dct["beta"] = lattice.beta
    dct["gamma"] = lattice.gamma

    # Flatten sites data
    dct["species"] = []
    dct["site_a"], dct["site_b"], dct["site_c"] = [], [], []
    dct["site_properties"] = []
    for site_dict in struct_dict["sites"]:
        dct["species"].append(site_dict["species"]["element"])
        dct["site_a"].append(site_dict["abc"][0])
        dct["site_b"].append(site_dict["abc"][1])
        dct["site_c"].append(site_dict["abc"][2])
        dct["site_properties"].append(site_dict["properties"])

    return dct


def _uncompress_struct_dict(dct: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """
    Reorganize structure data into object reconstruction-friendly organization in place.
    Exact reverse operation of `_compress_struct_dict()`.
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