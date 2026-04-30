"""Compute Elementary Metastability metric."""

import typing as tp
import itertools as itt
from collections import defaultdict, OrderedDict

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .metric_base import Metric, MetricsData
from src.utils import (
    GenMatPDEntry, GenMatStructure, PDEntryLike, StructureLike
)


class ElementMetastability(Metric):
    """
    Compute Element Metastability metric.
    
    Definition
    ----------
    Structures are considered metastable from elements if their formation energy
    is lower or equal to the convex hull energy between element references of their
    chemical space. The formation energy of a crystal is assumed to be the total energy
    of the unit cell divided by the number of atoms in it, and is measured in eV/atom.

    NOTE: Given structures must have their total energy attribute defined.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        ref_entries: list[PDEntry],
        workers: int | None = None,
        verbose: bool = False
    ) -> None:
        """
        Compute Element Metastability metric.

        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Metastability on.

        ref_unaries: list[PDEntry]
            Unary PDEntry objects of known formation energy to use as references.

        workers: int, optional
            Number of parallel processes to spawn for efficient computation. If not given,
            `tqdm.contrib.concurrent.process_map()` default is used. Pass 0 to disable
            `process_map()` and execute sequentially.

        verbose: bool
            Whether to print entries used to build convex hulls with their energies,
            and computed energies above hull from tested entries. Defaults to False.
        """
        super().__init__(structures, workers)
        assert all(struct.energy is not None for struct in self.structures)
        assert all(entry.composition.is_element for entry in ref_entries)
        self._check_lacking_elts(self.structures, ref_entries)

        self._sorted_structs = self.sort_by_chemical_system(self.structures)
        self.ref_entries = self.sort_by_chemical_system(ref_entries)
        self.verbose = verbose

        self._compute()

    @staticmethod
    def get_elements(structures: list[StructureLike] | list[PDEntryLike]) -> set[str]:
        """Get a set of all unique elements present in given structures or entries."""
        return set.union(*(set(struct.composition.get_el_amt_dict()) for struct in structures))

    def _check_lacking_elts(
        self,
        structures: list[GenMatStructure] | list[GenMatPDEntry],
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
    def sort_by_chemical_system(self, entries: list[GenMatPDEntry]) -> dict[str, list[GenMatPDEntry]]:
        ...
    @tp.overload
    def sort_by_chemical_system(self, entries: list[PDEntry]) -> dict[str, list[PDEntry]]:
        ...
    @tp.overload
    def sort_by_chemical_system(self, entries: list[GenMatStructure]) -> dict[str, list[GenMatStructure]]:
        ...
    def sort_by_chemical_system(self, entries: list) -> dict:
        """Sort entries by chemical system."""
        chemical_systems: dict[str, list[PDEntryLike]] = defaultdict(list[PDEntryLike])
        for entry in entries:
            chemical_systems[entry.composition.chemical_system].append(entry)
        return chemical_systems

    def is_stable(self, pd: PhaseDiagram, entry: GenMatStructure | PDEntryLike) -> bool:
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
        """
        if isinstance(entry, GenMatStructure):
            entry = entry.entry

        if self.verbose:
            pd_system = '-'.join(map(str, pd.elements))
            entry_name = entry.name or entry.formula
            print(f"Comparing entry {entry_name} to convex hull {pd_system}:")
            print(f"{entry.energy_per_atom=}")

        e_above_hull = pd.get_e_above_hull(entry, allow_negative=True, check_stable=False)
        is_stable = e_above_hull is not None and e_above_hull <= 0.0

        if self.verbose:
            print(f"{e_above_hull=}")
            print(f"element metastable: {is_stable}")

        return is_stable

    def get_stability(self, pd: PhaseDiagram, structure: GenMatStructure) -> MetricsData:
        """
        Compute ElementMetastability of a structure on a phase diagram.
        A computed `MetricsData` containing metric results is returned.
        Computed energy above hull is not saved.

        Parameters
        ----------
        pd: PhaseDiagram
            Reference diagram to use as convex hull.

        structure: GenMatStructure
            The structure to compute Stability on.

        Returns
        -------
        MetricsData
            Computed metric results.
        """
        return MetricsData(
            structure,
            is_element_metastable=self.is_stable(pd, structure),
            additional_data=self._get_metric_settings()
        )

    def _get_metric_settings(self) -> dict[str, tp.Any]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {}
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
    def elt_metastable_data(self) -> list[MetricsData]:
        """List of computed `MetricsData` containing element metastable structures."""
        return [data for data in self.computed_data if data.is_element_metastable]

    @property
    def elt_unstable_data(self) -> list[MetricsData]:
        """List of computed `MetricsData` containing element unstable structures."""
        return [data for data in self.computed_data if data.is_element_metastable is False]

    @property
    def elt_metastable_structs(self) -> list[GenMatStructure]:
        """List of element metastable structures."""
        return [data.structure for data in self.computed_data if data.is_element_metastable]

    @property
    def elt_unstable_structs(self) -> list[GenMatStructure]:
        """List of element unstable structures."""
        return [data.structure for data in self.computed_data if data.is_element_metastable is False]

    @property
    def elt_metastable_names(self) -> list[str]:
        """List of names of element metastable structures."""
        return [data.structure.name for data in self.computed_data if data.is_element_metastable]

    @property
    def elt_unstable_names(self) -> list[str]:
        """List of names of element unstable structures."""
        return [data.structure.name for data in self.computed_data if data.is_element_metastable is False]

    @property
    def elt_metastable_names_set(self) -> set[str]:
        """Set of names of element metastable structures."""
        return {data.structure.name for data in self.computed_data if data.is_element_metastable}

    @property
    def elt_unstable_names_set(self) -> set[str]:
        """Set of names of element unstable structures."""
        return {data.structure.name for data in self.computed_data if data.is_element_metastable is False}

    @property
    def elt_metastable_entries(self) -> list[GenMatPDEntry]:
        """List of `GenMatPDEntry` objects associated with element metastable structures."""
        return [data.structure.entry for data in self.computed_data if data.is_element_metastable]

    @property
    def elt_unstable_entries(self) -> list[GenMatPDEntry]:
        """List of `GenMatPDEntry` objects associated with element unstable structures."""
        return [data.structure.entry for data in self.computed_data if data.is_element_metastable is False]

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[GenMatStructure]] = OrderedDict(
            [
                ("Element metastable", self.elt_metastable_structs),
                ("Element unstable", self.elt_unstable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)