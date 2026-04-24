"""Compute Stability metric."""

import typing as tp
import typing_extensions as tpe
import itertools as itt
from collections import defaultdict, OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .metric_base import Metric, MetricsData, NoneMetric
from src.utils import (
    GenMatPDEntry, GenMatStructure, is_genmat_name, PDEntryLike, StructureLike,
    flatten, VisualIterator
)

# TODO: extend flexibility to allow for seemless computations of:
# - Element Metastability (Ehull <= 0.0, element refs only) => MetricsData.is_element_metastable
# - Metastability (Ehull <= stable_tol) => MetricsData.is_metastable
# - Stability (Ehull <= 0.0) => MetricsData.is_stable
class Stability(Metric):
    """
    Compute Stability metric.
    
    Definition
    ----------
    Structures are considered:
    - **Stable** if their relative formation energy is lower or equal to
    the convex hull energy of a reference dataset of known entries in their chemical space.
    - **Metastable** if the difference between their relative formation energy and the convex hull
    energy is lower than a threshold (generally 0.1 eV/atom).

    Notes
    -----
    - The formation energy of a crystal is assumed to be the total energy of the unit cell divided
    by the number of atoms in it, and is measured in eV/atom (See reference for theoretical
    discussion).
    - Given structures must have a defined energy value in the `energy` attribute.

    References
    ----------
    S. P. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from First Principles
    Calculations. Chem. Mater., 2008, 20(5), 1798-1807. doi:10.1021/cm702327g
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        ref_entries: list[PDEntry],
        stable_tol: float = 0.1,
        workers: int | None = None,
        verbose: bool = False
    ) -> None:
        """
        Compute Stability metric.

        Parameters
        ----------
        structures: list[GenMatStructure]
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
            and computed energies above hull from tested entries. Defaults to False.
        """
        super().__init__(structures, workers)
        if stable_tol < 0.0:
            raise ValueError(
                f"'stable_tol' must be positive or zero, got {stable_tol}."
            )
        if any(struct.energy is None for struct in self.structures):
            raise ValueError(
                "All structures must have a defined energy value "
                f"to compute {type(self).__name__}."
            )
        self._check_lacking_elts(structures, ref_entries)

        self.stable_tol = stable_tol
        self.verbose = verbose

        # Sort all entries by chemical system
        self._sorted_structs = self.sort_by_chemical_system(self.structures)
        self.ref_entries = self.sort_by_chemical_system(ref_entries)

        self._compute()

    @staticmethod
    def get_elements(structures: list[StructureLike] | list[PDEntry]) -> set[str]:
        """Get a set of all unique elements present in given structures or entries."""
        return set.union(*(set(struct.composition.get_el_amt_dict()) for struct in structures))

    def _check_lacking_elts(
        self,
        structures: list[StructureLike] | list[PDEntry],
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

    @tp.overload
    def sort_by_chemical_system(self, entries: list[GenMatStructure]) -> dict[str, list[GenMatStructure]]:
        ...
    @tp.overload
    def sort_by_chemical_system(self, entries: list[GenMatPDEntry]) -> dict[str, list[GenMatPDEntry]]:
        ...
    @tp.overload
    def sort_by_chemical_system(self, entries: list[PDEntry]) -> dict[str, list[PDEntry]]:
        ...
    def sort_by_chemical_system(self, entries: list) -> dict:
        """Sort entries by chemical system."""
        chemical_systems: dict[str, list[PDEntry]] = defaultdict(list[PDEntry])
        for entry in entries:
            chemical_systems[entry.composition.chemical_system].append(entry)
        return chemical_systems

    @staticmethod
    def _make_dummy_structure(entry: GenMatPDEntry) -> GenMatStructure:
        """Make a dummy GenMatStructure from a named PDEntry object."""
        species_iter = itt.chain.from_iterable(
            [[elt] * int(amt) for elt, amt in entry.composition.get_el_amt_dict().items()]
        )
        return GenMatStructure(
            name=entry.name,
            energy=entry.energy,
            lattice=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            species=list(species_iter),
            coords=[[0., 0., 0.]] * int(entry.composition.num_atoms)
        )

    @classmethod
    def from_entries(
        cls, entries: list[GenMatPDEntry], ref_entries: list[PDEntry], **kwargs
    ) -> tpe.Self:
        """
        Initialize from `GenMatPDEntry` objects instead of structures.
        This method wraps entries inside dummy structures to initialize the class.

        Parameters
        ----------
        entries: list[GenMatPDEntry]
            Entries to compute stability on.

        ref_entries: list[PDEntry]
            PDEntry objects of known formation energy to use as references.

        kwargs: Any
            Additional keyword arguments to pass to the constructor.

        Returns
        -------
        Stability
            Stability class instance.

        Warnings
        --------
        - The dummy structures are NOT valid in any physical sense. Original entries can be
        accessed after computation inside the dummy structures properties at the 'PDEntry'
        key, or using the `stable_entries` and `unstable_entries` attributes instead of
        `stable_structs` and `unstable_structs`.

        - The `write_result()` method with `verbose` set to `True` will show the list
        of the dummy structures and not the corresponding entries, making it unsuitable
        to use directly to print actual entries lists.
        """
        structures = [cls._make_dummy_structure(entry) for entry in entries]
        return cls(structures, ref_entries, **kwargs)

    def is_stable(
        self, pd: PhaseDiagram, entry: GenMatStructure | PDEntryLike, strict: bool = False
    ) -> bool:
        """
        Whether `entry` is stable compared to given convex hull.
        Compatible with standard `PDEntry` and `GenMatPDEntry` objects, but only return
        a bool instead of a computed `MetricsData` object with metric results.
        Computed energy above hull is not saved.

        Parameters
        ----------
        pd: PhaseDiagram
            Reference diagram to use as convex hull.

        entry: GenMatStructure | PDEntry
            The entry to compute stability on.

        strict: bool
            Whether to compute strict stability or metastability using initialized tolerance.
            Defaults to False.
        """
        if isinstance(entry, GenMatStructure):
            entry = entry.entry

        if self.verbose:
            pd_system = '-'.join(map(str, pd.elements))
            entry_name = entry.name or entry.formula
            print(f"Comparing entry {entry_name} to convex hull {pd_system}:")
            print(f"{entry.energy_per_atom=}")

        e_above_hull = pd.get_e_above_hull(entry, allow_negative=True, check_stable=False)
        stable_tol = 0.0 if strict else self.stable_tol
        is_stable = e_above_hull is not None and e_above_hull <= stable_tol

        if self.verbose:
            stable_type = "stable" if strict else "metastable"
            print(f"{e_above_hull=}")
            print(f"{stable_type}: {is_stable}")

        return is_stable

    def get_stability(self, pd: PhaseDiagram, structure: GenMatStructure) -> MetricsData:
        """
        Compute Stability of a structure on a phase diagram, within initialized tolerance.
        The structure `energy_above_hull` attribute is populated with computed energy
        above hull, and a computed `MetricsData` containing metric results is returned.

        Parameters
        ----------
        pd: PhaseDiagram
            Reference diagram to use as convex hull.

        structure: GenMatStructure
            The structure to compute Stability on.

        Returns
        -------
        MetricsData
            Computed metric results containing the structure with updated energy above hull.
        """
        structure.compute_energy_above_hull(pd, allow_negative=True, check_stable=False)
        is_metastable = (
            structure.energy_above_hull is not None and
            structure.energy_above_hull <= self.stable_tol
        )
        is_stable = (
            structure.energy_above_hull is not None and
            structure.energy_above_hull <= 0.0
        )
        return MetricsData(
            structure,
            is_metastable=is_metastable,
            is_stable=is_stable,
            additional_data=self._get_metric_settings()
        )

    def _get_metric_settings(self) -> dict[str, tp.Any]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "stable_tol": self.stable_tol
            }
        }

    def _compute_system(
        self, system: str, structs: list[GenMatStructure]
    ) -> list[MetricsData]:
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
        return [self.get_stability(pd, struct) for struct in structs]

    def _compute(self) -> None:
        if self.workers == 0:
            self._computed_data: list[MetricsData] = list(itt.chain.from_iterable(
                self._sequential_compute(
                    self._compute_system, self._sorted_structs.items(), unpack=True
                )
            ))
        else:
            self._computed_data: list[MetricsData] = list(itt.chain.from_iterable(
                self._parallel_compute(
                    self._compute_system, zip(*self._sorted_structs.items()), unpack=True
                )
            ))

    @property
    def computed_data(self) -> list[MetricsData]:
        """List of all computed `MetricsData` objects."""
        return self._computed_data

    @property
    def stable_data(self) -> list[MetricsData]:
        """List of computed `MetricsData` containing a strictly stable structure."""
        return [data for data in self.computed_data if data.is_stable]

    @property
    def metastable_data(self) -> list[MetricsData]:
        """List of computed `MetricsData` containing a metastable structure."""
        return [data for data in self.computed_data if data.is_metastable]

    @property
    def unstable_data(self) -> list[MetricsData]:
        """List of computed `MetricsData` containing an unstable structure."""
        return [data for data in self.computed_data if data.is_metastable is False]

    @property
    def stable_structs(self) -> list[GenMatStructure]:
        """List of strictly stable structures."""
        return [data.structure for data in self.computed_data if data.is_stable]

    @property
    def metastable_structs(self) -> list[GenMatStructure]:
        """List of metastable structures."""
        return [data.structure for data in self.computed_data if data.is_metastable]

    @property
    def unstable_structs(self) -> list[GenMatStructure]:
        """List of unstable structures."""
        return [data.structure for data in self.computed_data if data.is_metastable is False]

    @property
    def stable_names(self) -> list[str]:
        """List of names of strictly stable structures."""
        return [data.structure.name for data in self.computed_data if data.is_stable]

    @property
    def metastable_names(self) -> list[str]:
        """List of names of metastable structures."""
        return [data.structure.name for data in self.computed_data if data.is_metastable]

    @property
    def unstable_names(self) -> list[str]:
        """List of names of unstable structures."""
        return [data.structure.name for data in self.computed_data if data.is_metastable is False]

    @property
    def stable_names_set(self) -> set[str]:
        """Set of names of strictly stable structures."""
        return {data.structure.name for data in self.computed_data if data.is_stable}

    @property
    def metastable_names_set(self) -> set[str]:
        """Set of names of metastable structures."""
        return {data.structure.name for data in self.computed_data if data.is_metastable}

    @property
    def unstable_names_set(self) -> set[str]:
        """Set of names of unstable structures."""
        return {data.structure.name for data in self.computed_data if data.is_metastable is False}

    @property
    def stable_entries(self) -> list[GenMatPDEntry]:
        """List of the `GenMatPDEntry` associated with stable structures."""
        return [data.structure.entry for data in self.computed_data if data.is_stable]

    @property
    def metastable_entries(self) -> list[GenMatPDEntry]:
        """List of the `GenMatPDEntry` associated with metastable structures."""
        return [data.structure.entry for data in self.computed_data if data.is_metastable]

    @property
    def unstable_entries(self) -> list[GenMatPDEntry]:
        """List of the `GenMatPDEntry` associated with unstable structures."""
        return [data.structure.entry for data in self.computed_data if data.is_metastable is False]

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[GenMatStructure]] = OrderedDict(
            [
                ("Stable", self.stable_structs),
                ("Metastable", self.metastable_structs),
                ("Unstable", self.unstable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
