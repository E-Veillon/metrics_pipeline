"""GenMat pipeline specific functionnalities."""

import re
import typing as tp
import typing_extensions as tpe
from pathlib import Path

from pymatgen.core.structure import Structure, Composition
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .spg_data import Spacegroup
from .visual_iterator import VisualIterator
from .common_asserts import check_type, check_num_value
from src.computations.local import get_density
from src.io import JsonLoader, JsonWriter, PathLike


_ELT_WITH_INDEX_REGEX = r"[A-Z][a-z]?\d*"
_GENMAT_NAME_REGEX = fr"^(\d+)_(?:{_ELT_WITH_INDEX_REGEX}|\((?:{_ELT_WITH_INDEX_REGEX})+\)\d*)+$"

def is_genmat_name(name: str) -> bool:
    """
    Whether given name follows GenMat naming convention.
    """
    return isinstance(name, str) and re.fullmatch(_GENMAT_NAME_REGEX, name) is not None


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
    match = re.match(_GENMAT_NAME_REGEX, name)
    assert match is not None
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
    A small extension of pymatgen's `Structure` objects to store GenMat name,
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

    energy_above_hull: float | None
        The energy above hull of the structure, in eV/atom.
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
        A small extension of pymatgen's `Structure` objects to store GenMat name,
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
        match = re.match(_GENMAT_NAME_REGEX, self.name)
        assert match is not None
        return int(match.group(1))

    @property
    def entry(self) -> GenMatPDEntry:
        """The `GenMatPDEntry` object associated with the structure."""
        if self.energy is None:
            raise ValueError(f"Energy not defined for structure {self.name}.")
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

    def compute_energy_above_hull(self, pdiagram: PhaseDiagram, **kwargs) -> None:
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

    def as_dict(self, *args, **kwargs) -> dict:
        """Get a dictionary representation of the GenMatStructure object."""
        struct_dct = super().as_dict(*args, **kwargs)
        del struct_dct["@module"]
        del struct_dct["@class"]
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

    @classmethod
    def from_dict(cls, dct: dict) -> tpe.Self:
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
    structures: list[Structure], special_keys: dict[str, tp.Any] | None = None
) -> list[GenMatStructure]:
    """
    Generate GenMatStructure objects from a list of Structure objects.
    The GenMat name is generated based on the position of each structure in the list,
    and stored in the 'name' attribute of the GenMatStructure object.
    """
    genmat_structures = [
        GenMatStructure.from_structure(
            name=f"{idx}_{structure.reduced_formula}",
            structure=structure,
            special_keys=special_keys or {}
        ) for idx, structure in enumerate(structures)
    ]
    return genmat_structures


class GenMatFile:
    """
    Lightweight wrapper to load and write JSON files containing structure data for
    the GenMat library.
    """
    _genmat_keys = (
        "name", "special_keys", "spacegroup", "energy", "energy_above_hull", "structure"
    )

    def __init__(self, filepath: PathLike | None = None) -> None:
        """
        Lightweight wrapper to load and write JSON files containing structure data for
        the GenMat library.

        Parameters
        ----------
        filepath: str | Path, optional
            The JSON file to read data from. If not given, initialize an empty file object
            that can be filled with data to write.
        """
        if filepath is None:
            self._data = {}
            for key in self._genmat_keys:
                self._data[key] = []
        else:
            data = JsonLoader(filepath).load_as_dict()
            assert data.get("@module") == type(self).__module__
            assert data.get("@class") == type(self).__name__
            assert all(key in data for key in self._genmat_keys)
            assert len({len(data[key]) for key in self._genmat_keys}) == 1

            self._data = {k: data[k] for k in self._genmat_keys}

    def __getattr__(self, name: str) -> tp.Any:
        if (attr_data:=self._data.get(name)) is not None:
            return attr_data
        raise AttributeError(f"{type(self).__name__!r} has no attribute {name!r}.")

    def parse_structures(self) -> list[GenMatStructure]:
        """Build GenMatStructure objects from file data."""
        return [
            GenMatStructure.from_dict(dict(zip(self._genmat_keys, values)))
            for values in zip(*(self._data[k] for k in self._genmat_keys))
        ]

    def add_structure(self, structure: GenMatStructure) -> None:
        """Add a GenMatStructure data to the file."""
        dct = structure.as_dict()
        for key in self._genmat_keys:
            self._data[key].append(dct[key])

    def add_structures(self, structures: list[GenMatStructure]) -> None:
        """Add several GenMatStructures to the file. Order is preserved."""
        for key in self._genmat_keys:
            self._data[key].extend(getattr(struct, key) for struct in structures)

    def write_file(self, filepath: PathLike) -> None:
        """Write stored data to a compact JSON file."""
        data = self._data.copy()
        data["@module"] = type(self).__module__
        data["@class"] = type(self).__name__
        JsonWriter(filepath, data).write_as_dict()
