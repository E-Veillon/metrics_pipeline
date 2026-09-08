"""GenMat pipeline specific functionalities."""

import re
import typing as tp
import typing_extensions as tpe

from pymatgen.core.structure import Structure, Composition
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .spg_data import Spacegroup
from .common_asserts import check_type


_ELT_WITH_INDEX_REGEX = r"[A-Z][a-z]?\d*"


class GenMatName(str):
    """
    A string subclass to represent GenMat names.
    It is validated at initialization to follow GenMat naming convention.
    """
    _GENMAT_NAME_PATTERN = re.compile(
        fr"^(\d+)_({_ELT_WITH_INDEX_REGEX}|\((?:{_ELT_WITH_INDEX_REGEX})+\)\d*)+$"
    )

    def __init__(self, name: str) -> None:
        """"""
        check_type(name, "name", (str,))
        match = self._GENMAT_NAME_PATTERN.fullmatch(name)
        if match is None:
            raise ValueError(
                f"{name!r} does not follow GenMat naming convention: '[index]_[formula]'."
            )
        self._name = name
        self._index = int(match.group(1))
        self._formula = match.group(2)

    @classmethod
    def is_genmat_name(cls, name: str) -> bool:
        """Whether given string follows GenMat naming convention."""
        return isinstance(name, str) and cls._GENMAT_NAME_PATTERN.fullmatch(name) is not None

    @property
    def name(self) -> str:
        """The full GenMat name as a string."""
        return self._name

    @property
    def index(self) -> int:
        """The index part of the GenMat name."""
        return self._index

    @property
    def formula(self) -> str:
        """The formula part of the GenMat name."""
        return self._formula


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

    name: GenMatName
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
        name: GenMatName,
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

        name: GenMatName
            A name for the entry. This is the string shown in the phase diagrams.
            Must follow GenMat naming convention.

        energy_above_hull: float, optional
            Computed energy above hull of the entry when compared to a `PhaseDiagram` object.

        attribute: MSONable, optional
            A arbitrary attribute. Can be used to specify that the
            entry is a newly found compound, or to specify a particular label for
            the entry, etc. An attribute can be anything but must be MSONable.
        """
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

    def as_pdentry(self) -> PDEntry:
        """Rebuild standard `PDEntry` object from data."""
        return PDEntry(self.composition, self.energy, self.name, self.attribute)
    
    @classmethod
    def from_pdentry(
        cls, entry: PDEntry, name: GenMatName | None = None, energy_above_hull: float | None = None
    ) -> tpe.Self:
        """
        Build a `GenMatPDEntry` object from a standard `PDEntry` object.
        
        Parameters
        ----------
        name: GenMatName, optional
            The GenMat name of the entry. If not given, current name of the entry
            will be checked and kept if it follows GenMat naming convention.

        energy_above_hull: float, optional
            Computed energy above hull of the entry when compared to a `PhaseDiagram` object.
        """
        return cls(
            composition=entry.composition,
            energy=entry.energy,
            name=name if name is not None else GenMatName(entry.name),
            energy_above_hull=energy_above_hull,
            attribute=entry.attribute
        )


class GenMatStructure(Structure):
    """
    An extension of pymatgen's `Structure` objects to store GenMat name,
    optional saved CIF labels, and useful properties for metrics computations.
    Only extension attributes are shown below, but the class also inherits all
    `Structure` attributes.

    Attributes
    ----------
    name: GenMatName
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
        name: GenMatName,
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
        name: GenMatName
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
        return self.name.index

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
        return sum(s.atomic_mass.to("g") for s in self.species) / self.volume*1e-24

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
            "energy_above_hull": self.energy_above_hull
        } | struct_dct
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
        dct_copy = dct.copy()
        dct_copy.pop("@module", None)
        dct_copy.pop("@class", None)

        # Build base structure
        struct_dict = {
            "@module": Structure.__module__,
            "@class": Structure.__name__,
            "lattice": dct_copy.pop("lattice"),
            "sites": dct_copy.pop("sites"),
            "charge": dct_copy.pop("charge", None),
            "properties": dct_copy.pop("properties", None)
        }
        structure = Structure.from_dict(**struct_dict)

        # Build spacegroup if defined
        spg = None
        if "spacegroup" in dct_copy:
            spg = Spacegroup.from_dict(dct_copy.pop("spacegroup"))

        return cls.from_structure(structure=structure, spacegroup=spg, **dct_copy)

    @classmethod
    def from_structure(cls, name: GenMatName, structure: Structure, **kwargs) -> tpe.Self:
        """
        Build a `GenMatStructure` from an existing `Structure` object.
        This conversion is only meant as an extension of the original structure to
        the features of `GenMatStructure` class. The `Structure` data will be kept as-is
        without any additional check or transformation purposefully for data fidelity.

        Parameters
        ----------
        name: GenMatName
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
            name=GenMatName(f"{idx}_{structure.reduced_formula}"),
            structure=structure,
            special_keys={k: structure.properties.get(k) for k in special_keys}
        ) for idx, structure in enumerate(structures)
    ]
    return genmat_structures
