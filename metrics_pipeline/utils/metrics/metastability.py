"""Compute Elementary Metastability metric."""

import itertools as itt
from collections import defaultdict

from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

from .metric_base import Metric


class ElementaryMetastability(Metric):
    """
    Compute Elementary Metastability metric.
    
    Definition
    ----------
    Structures are considered metastable from elements if their relative formation energy
    is lower or equal to the convex hull energy between each unary crystal of their
    chemical space. The formation energy of a crystal is assumed to be the total energy
    of the unit cell divided by the number of atoms in it, and is measured in eV/atom.

    NOTE: Given structures (both computed and references) must have their total energy
    in eV stored in their properties under the 'energy' key.
    """
    def __init__(self, structures: list[Structure], ref_unaries: list[Structure]) -> None:
        """
        Compute Elementary Metastability metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Metastability on.

        ref_unaries: list[Structure]
            Unary structures of known formation energy to use as elementary references.
        """
        super().__init__(structures)

        assert all(struct.properties.get("energy") is not None for struct in self.structures)
        assert all(struct.properties.get("energy") is not None for struct in ref_unaries)
        assert all(len(struct.composition) == 1 for struct in ref_unaries)
        used_elts = self.get_elements(self.structures)
        ref_elts = self.get_elements(ref_unaries)
        lacking_elts = list(used_elts - ref_elts)
        assert not lacking_elts, ValueError(
            "Following elements are present in structures but lacking in references: "
            f"{', '.join(sorted(lacking_elts))}."
        )

        self.ref_unaries = ref_unaries

        self._compute()

    @staticmethod
    def get_elements(structures: list[Structure]) -> set[str]:
        """Get a set of all unique elements present in given structures."""
        return set( # eliminate duplicate elements
            itt.chain.from_iterable( # flatten from list[list[str]] to list[str]
                list(map(str, struct.composition.elements)) for struct in structures
            )
        )

    @staticmethod
    def is_stable(pd: PhaseDiagram, entry: PDEntry) -> bool:
        """Whether entry is stable compared to the convex hull."""
        e_above_hull = pd.get_e_above_hull(entry, allow_negative=True, check_stable=False)
        if e_above_hull is None:
            return False
        return e_above_hull <= 0.0

    def _compute(self) -> None:
        self._metastable_structs: list[Structure] = []
        self._unstable_structs: list[Structure] = []
        chemical_systems = defaultdict(list[PDEntry])
        # Group entries to compute into distinct chemical systems
        for struct in self.structures:
            chemical_systems[struct.chemical_system].append(
                PDEntry(struct.composition, struct.properties["energy"], attribute=struct)
            )
        # Iterate through each chemical system to compute metastable and unstable
        # structures in this system
        for system, entries in chemical_systems.items():
            ref_entries = [
                PDEntry(struct.composition, struct.properties["energy"])
                for struct in self.ref_unaries if struct.elements[0].symbol in system
            ]
            pd = PhaseDiagram(ref_entries)
            self._metastable_structs.extend(
                [entry.attribute for entry in entries if self.is_stable(pd, entry)] # type: ignore
            )
            self._unstable_structs.extend(
                [entry.attribute for entry in entries if not self.is_stable(pd, entry)] # type: ignore
            )

    @property
    def metastable_structs(self) -> list[Structure]:
        return self._metastable_structs

    @property
    def unstable_structs(self) -> list[Structure]:
        return self._unstable_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        return self._write_filter_metric_result(
            self._metastable_structs, self._unstable_structs, filename, verbose
        )