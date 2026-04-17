"""Compute Stability metric."""

import typing as tp
import typing_extensions as tpe
import itertools as itt
from collections import defaultdict, OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .metric_base import Metric
from src.utils import flatten, VisualIterator

class Stability(Metric):
    """
    Compute Stability metric.
    
    Definition
    ----------
    Structures are considered stable if their relative formation energy is lower or equal to
    the convex hull energy of a reference dataset of known entries in their chemical space.
    The formation energy of a crystal is assumed to be the total energy of the unit cell divided
    by the number of atoms in it, and is measured in eV/atom.

    NOTE: Given structures (both computed and references) must have their total energy
    in eV stored in their properties under the 'energy' key.
    """
    _struct_attr = f"{__qualname__}_structure"
    _delta_e_attr = f"{__qualname__}_e_above_hull"
    def __init__(
        self,
        structures: list[Structure],
        ref_entries: list[PDEntry],
        stable_tol: float = 0.1,
        workers: int | None = None,
        verbose: bool = False
    ) -> None:
        """
        Compute Stability metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Stability on.

        ref_entries: list[PDEntry]
            PDEntry objects of known formation energy to use as references.

        stable_tol: float
            Tolerance of energy above the hull to consider the structure as stable.
            Defaults to 0.1 eV/atom.

        workers: int, optional
            Number of parallel processes to spawn for efficient computation. If not given,
            `tqdm.contrib.concurrent.process_map()` default is used. Pass 0 to disable
            `process_map()` and execute sequentially.

        verbose: bool
            Whether to print entries used to build convex hulls with their energies,
            and computed energies above hull from tested entries.
        """
        super().__init__(structures)
        self._init_settings(stable_tol, workers, verbose)
        self._check_lacking_elts(structures, ref_entries)
        if not all(
            isinstance(struct.properties.get("header"), str)
            for struct in self.structures
        ):
            raise ValueError(
                "All structures must have a 'header' string property "
                "to be used as an identifier in results."
            )
        if not all(
            isinstance(struct.properties.get("energy"), float)
            for struct in self.structures
        ):
            raise ValueError(
                "All structures must have an 'energy' value in properties."
            )
        # Convert all structures to lightweight PDEntry objects
        entries = [
            PDEntry(
                struct.composition, struct.properties["energy"], struct.properties["header"]
            ) for struct in self.structures
        ]
        # Sort all entries by chemical system
        self.entries = self.sort_by_chemical_system(entries)
        self.ref_entries = self.sort_by_chemical_system(ref_entries)

        self._compute()

    @staticmethod
    def _make_dummy_structure(entry: PDEntry) -> Structure:
        """Make a dummy structure to store a PDEntry object in its properties."""
        species_iter = itt.chain.from_iterable(
            [[elt] * int(amt) for elt, amt in entry.composition.get_el_amt_dict().items()]
        )
        return Structure(
            lattice=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            species=list(species_iter),
            coords=[[0., 0., 0.]] * int(entry.composition.num_atoms),
            properties={"header": entry.name, "energy": entry.energy}
        )

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
            PDEntry objects of known formation energy to use as references.

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
        if not all(isinstance(entry.name, str) for entry in entries):
            raise ValueError(
                "All entries must have a defined string name to be used "
                "as an identifier in results."
            )
        structures = [cls._make_dummy_structure(entry) for entry in entries]
        return cls(structures, ref_entries, **kwargs)

    @staticmethod
    def get_elements(structures: list[Structure] | list[PDEntry]) -> set[str]:
        """Get a set of all unique elements present in given structures or entries."""
        return set.union(*(set(struct.composition.get_el_amt_dict()) for struct in structures))

    def _init_settings(self, stable_tol: float, workers: int | None, verbose: bool) -> None:
        """Verify and initialize settings arguments."""
        if stable_tol < 0.0:
            raise ValueError(
                f"'stable_tol' must be positive or zero, got {stable_tol}."
            )
        if workers is not None and workers < 0:
            raise ValueError(
                f"'workers must be positive or zero, got {workers}."
            )
        self.stable_tol = stable_tol
        self.workers = workers
        self.verbose = verbose

    def _check_lacking_elts(
        self,
        structures: list[Structure] | list[PDEntry],
        ref_entries: list[PDEntry]
    ) -> None:
        """Verify that no element present in candidates is lacking in references."""
        used_elts = self.get_elements(structures)
        ref_elts = self.get_elements(ref_entries)
        lacking_elts = used_elts - ref_elts
        if lacking_elts:
            raise ValueError(
                "Following elements are present in structures but lacking in references: "
                f"{', '.join(sorted(lacking_elts))}."
        )

    def sort_by_chemical_system(self, entries: list[PDEntry]) -> dict[str, list[PDEntry]]:
        """Sort entries by chemical system."""
        chemical_systems: dict[str, list[PDEntry]] = defaultdict(list[PDEntry])
        for entry in entries:
            chemical_systems[entry.composition.chemical_system].append(entry)
        return chemical_systems

    def is_stable(self, pd: PhaseDiagram, entry: PDEntry) -> bool:
        """Whether entry is stable compared to the convex hull, within initialized tolerance."""
        if self.verbose:
            pd_system = '-'.join(map(str, pd.elements))
            print(f"Comparing entry {entry} to convex hull {pd_system}:")
            print(f"{entry.energy_per_atom=}")
        e_above_hull = pd.get_e_above_hull(entry, allow_negative=True, check_stable=False)
        if entry.attribute is None:
            entry.attribute = {}
        entry.attribute[self._delta_e_attr] = e_above_hull # type: ignore
        if self.verbose:
            print(f"{e_above_hull=}")
        if e_above_hull is None or e_above_hull > self.stable_tol:
            if self.verbose:
                print(f"stable: False")
            return False
        if self.verbose:
            print("stable: True")
        return True

    def _compute_system(
        self, system: str, entries: list[PDEntry]
    ) -> tuple[list[PDEntry], list[PDEntry]]:
        """Compute stability inside a chemical system."""
        if self.verbose:
            print(f"Computing system {system}")
        # Get relevant reference entries for this chemical system
        system_set = set(system.split("-"))
        ref_entries = list(
            itt.chain.from_iterable(
                self.ref_entries[ref_system] for ref_system in self.ref_entries
                if set(ref_system.split("-")).issubset(system_set)
            )
        )
        # Build the reference phase diagram
        pd = PhaseDiagram(ref_entries)
        if self.verbose:
            print("Built phase diagram with following stable entries:")
            hull_entries = tp.cast(tuple[PDEntry, ...], pd.qhull_entries)
            for entry in hull_entries:
                print(f"{entry=}, {entry.energy_per_atom=} eV/atom")
        # Compute energy above hull of evaluated entries
        stable_entries, unstable_entries = [], []
        for entry in entries:
            (stable_entries if self.is_stable(pd, entry) else unstable_entries).append(entry)

        return stable_entries, unstable_entries

    def _compute(self) -> None:
        desc="Computing Stability"
        if self.workers == 0:
            self._stable_entries: list[PDEntry] = []
            self._unstable_entries: list[PDEntry] = []
            iterator = VisualIterator(
                self.entries.items(), desc=desc, unit="computed", percent=True
            )
            for system, entries in iterator:
                stable_entries, unstable_entries = self._compute_system(system, entries)
                self._stable_entries.extend(stable_entries)
                self._unstable_entries.extend(unstable_entries)
        else:
            stable_entries, unstable_entries = zip(
                *process_map(
                    self._compute_system,
                    *zip(*self.entries.items()),
                    max_workers=self.workers,
                    chunksize=min(10, len(self.entries) // 100 + 1),
                    desc=desc
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
    def stable_names(self) -> list[str]:
        """List of identifiers of stable structures."""
        return [entry.name for entry in self.stable_entries]

    @property
    def unstable_names(self) -> list[str]:
        """List of identifiers of unstable structures."""
        return [entry.name for entry in self.unstable_entries]

    @property
    def stable_names_set(self) -> set[str]:
        """Set of identifiers of stable structures."""
        return set(entry.name for entry in self.stable_entries)

    @property
    def unstable_names_set(self) -> set[str]:
        """Set of identifiers of unstable structures."""
        return set(entry.name for entry in self.unstable_entries)

    @property
    def stable_structs(self) -> list[Structure]:
        """List of stable structures."""
        return [
            struct for struct in self.structures
            if struct.properties["header"] in self.stable_names_set
        ]

    @property
    def unstable_structs(self) -> list[Structure]:
        """List of unstable structures."""
        return [
            struct for struct in self.structures
            if struct.properties["header"] in self.unstable_names_set
        ]

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Stable", self.stable_structs),
                ("Unstable", self.unstable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
