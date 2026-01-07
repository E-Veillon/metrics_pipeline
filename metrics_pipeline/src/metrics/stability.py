"""Compute Stability metric."""

import typing as tp
import typing_extensions as tpe
import itertools as itt
from collections import defaultdict, OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .metric_base import Metric
from src.utils import flatten

class Stability(Metric):
    """
    Compute Stability metric.
    
    Definition
    ----------
    Structures are considered stable if their relative formation energy is lower or equal to
    the convex hull energy of a reference dataset of known structures in their chemical space.
    The formation energy of a crystal is assumed to be the total energy of the unit cell divided
    by the number of atoms in it, and is measured in eV/atom.

    NOTE: Given structures (both computed and references) must have their total energy
    in eV stored in their properties under the 'energy' key.
    """
    _struct_attr = f"{__name__}_structure"
    _delta_e_attr = f"{__name__}_e_above_hull"
    def __init__(
        self,
        structures: list[Structure],
        ref_structs: list[Structure],
        stable_tol: float = 0.1,
        workers: int | None = None
    ) -> None:
        """
        Compute Stability metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Stability on.

        ref_structs: list[Structure]
            Structures of known formation energy to use as references.

        stable_tol: float
            Tolerance of energy above the hull to consider the structure as stable.
            Defaults to 0.1 eV/atom.

        workers: int, optional
            Number of parallel processes to spawn for efficient computation. If not given,
            `tqdm.contrib.concurrent.process_map()` default is used. Pass 0 to disable
            `process_map()` and execute sequentially.
        """
        super().__init__(structures)

        assert all(struct.properties.get("energy") is not None for struct in self.structures)
        assert all(struct.properties.get("energy") is not None for struct in ref_structs)
        used_elts = self.get_elements(self.structures)
        ref_elts = self.get_elements(ref_structs)
        lacking_elts = list(used_elts - ref_elts)
        assert not lacking_elts, ValueError(
            "Following elements are present in structures but lacking in references: "
            f"{', '.join(sorted(lacking_elts))}."
        )
        assert stable_tol >= 0.0, ValueError(
            f"'stable_tol' must be positive or zero, got {stable_tol}."
        )
        if workers is not None:
            assert workers >= 0, ValueError(
                f"'workers must be positive or zero, got {workers}."
            )
        self.ref_structs = ref_structs
        self.stable_tol = stable_tol
        self.workers = workers

        chemical_systems: dict[str, list[PDEntry]] = defaultdict(list[PDEntry])
        # Group entries to compute into distinct chemical systems
        for struct in self.structures:
            internal_entry = struct.properties.get("PDEntry")
            if internal_entry is None:
                entry = PDEntry(struct.composition, struct.properties["energy"])
                entry.attribute = {
                    "orig_attribute": entry.attribute,
                    self._struct_attr: struct
                }
                chemical_systems[struct.chemical_system].append(entry)
            else:
                assert isinstance(internal_entry, PDEntry), RuntimeError(
                    "Type checker assertion."
                )
                internal_entry.attribute = {
                    "orig_attribute": internal_entry.attribute,
                    self._struct_attr: struct
                }
                chemical_systems[struct.chemical_system].append(internal_entry)

        self.entries_dict = chemical_systems
        self.ref_entries = [
            PDEntry(struct.composition, struct.properties["energy"])
            if struct.properties.get("PDEntry") is None else struct.properties["PDEntry"]
            for struct in self.ref_structs
        ]

        self._compute()

        self._stable_structs: list[Structure] = [
            entry.attribute[self._struct_attr] # type: ignore
            for entry in self._stable_entries
        ]
        self._unstable_structs: list[Structure] = [
            entry.attribute[self._struct_attr] # type: ignore
            for entry in self._stable_entries
        ]

    @classmethod
    def from_entries(
        cls, entries: list[PDEntry], ref_entries: list[PDEntry], **kwargs
    ) -> tpe.Self:
        """
        Initialize from PDEntry objects instead of structures.
        This method wraps entries inside dummy structures to initialize the class.

        Parameters
        ----------
        entries: list[PDEntry]
            Entries to compute stability on.

        ref_entries: list[PDEntry]
            Entries of known formation energy to use as references.

        kwargs: Any
            Additional keyword arguments to pass to the constructor.

        Returns
        -------
        Stability
            Stability class instance.

        Warning
        -------
        - The dummy structures are NOT valid in any physical sense. Original entries can be
        accessed after computation inside the dummy structures properties at the 'PDEntry'
        key, or using the `stable_entries` and `unstable_entries` attributes instead of
        `stable_structs` and `unstable_structs`.

        - The `write_result()` method with `verbose` set to `True` will show the list
        of the dummy structures and not the corresponding entries, making it unsuitable
        to use directly to print actual entries lists.
        """
        structures = [
            Structure(
                lattice=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                species=flatten([[elt] * int(nbr) for elt, nbr in entry.composition.items()]),
                coords=[[0., 0., 0.] for _ in range(int(entry.composition.num_atoms))],
                properties={"energy": entry.energy, "PDEntry": entry}
            ) for entry in entries
        ]
        ref_structs = [
            Structure(
                lattice=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                species=flatten([[elt] * int(nbr) for elt, nbr in entry.composition.items()]),
                coords=[[0., 0., 0.] for _ in range(int(entry.composition.num_atoms))],
                properties={"energy": entry.energy, "PDEntry": entry}
            ) for entry in ref_entries
        ]
        stability = cls(structures, ref_structs, **kwargs)
        return stability

    @staticmethod
    def get_elements(structures: list[Structure]) -> set[str]:
        """Get a set of all unique elements present in given structures."""
        return set( # eliminate duplicate elements
            itt.chain.from_iterable( # flatten from list[list[str]] to list[str]
                list(map(str, struct.composition.elements)) for struct in structures
            )
        )

    def is_stable(self, pd: PhaseDiagram, entry: PDEntry) -> bool:
        """Whether entry is stable compared to the convex hull, within initialized tolerance."""
        e_above_hull = pd.get_e_above_hull(entry, allow_negative=True, check_stable=False)
        entry.attribute[self._delta_e_attr] = e_above_hull # type: ignore
        if e_above_hull is None:
            return False
        return e_above_hull <= self.stable_tol

    def _compute_system(
        self, system: str, entries: list[PDEntry]
    ) -> tuple[list[PDEntry], list[PDEntry]]:
        """Compute stability inside a chemical system."""
        system_set = set(system.split("-"))
        ref_entries = [
            entry for entry in self.ref_entries
            if entry.composition.chemical_system_set.issubset(system_set)
        ]
        pd = PhaseDiagram(ref_entries)
        stable_structs = [entry for entry in entries if self.is_stable(pd, entry)]
        unstable_structs = [entry for entry in entries if not self.is_stable(pd, entry)]

        return stable_structs, unstable_structs

    def _compute(self) -> None:
        if self.workers == 0:
            self._stable_entries: list[PDEntry] = []
            self._unstable_entries: list[PDEntry] = []
            for system, entries in self.entries_dict.items():
                stable_entries, unstable_entries = self._compute_system(system, entries)
                self._stable_entries.extend(stable_entries)
                self._unstable_entries.extend(unstable_entries)
        else:
            stable_entries, unstable_entries = zip(
                *process_map(
                    self._compute_system,
                    *zip(*self.entries_dict.items()),
                    max_workers=self.workers,
                    chunksize=min(10, len(self.entries_dict) // 100 + 1),
                    desc="Computing Stability"
                )
            )
            self._stable_entries = flatten(stable_entries)
            self._unstable_entries = flatten(unstable_entries)

    @property
    def stable_entries(self) -> list[PDEntry]:
        """List of stable phase diagram entries."""
        return self._stable_entries

    @property
    def unstable_entries(self) -> list[PDEntry]:
        """List of unstable phase diagram entries."""
        return self._unstable_entries

    @property
    def stable_structs(self) -> list[Structure]:
        """List of stable structures."""
        return self._stable_structs

    @property
    def unstable_structs(self) -> list[Structure]:
        """List of unstable structures."""
        return self._unstable_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Stable", self._stable_structs),
                ("Unstable", self._unstable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)