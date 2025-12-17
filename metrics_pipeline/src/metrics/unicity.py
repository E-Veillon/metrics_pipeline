"""Compute Unicity metric."""

import itertools as itt
from collections import OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from .metric_base import Metric
from .hash_matcher import group_compositions

class Unicity(Metric):
    """
    Compute Unicity metric.
    
    Definition
    ----------
    Structures are uniques if they are not an equivalent representation of a previous structure
    in given list.
    """
    def __init__(
        self,
        structures: list[Structure],
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5.0,
        workers: int | None = None
    ) -> None:
        """
        Compute Unicity metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Unicity on.

        ltol: float
            Fractional length tolerance for structures matching. Default is 0.2.

        stol: float
            Site tolerance for structures matching. Defined as the fraction of the
            average free length per atom, i.e. ( V / Nsites ) ** (1/3). Default is 0.3.

        angle_tol: float
            Angle tolerance for structures matching in degrees. Default is 5.0 degrees.

        workers: int, optional
            Number of parallel processes to spawn for high throughput structure matching.
            If not given, will use default of `tqdm.contrib.concurrent.process_map()`.
            pass 0 to disable `process_map()` and do it sequentially.
        """
        super().__init__(structures)
        self.matcher = StructureMatcher(
            ltol=ltol, stol=stol, angle_tol=angle_tol, scale=False, attempt_supercell=True
        )
        if workers is not None:
            assert isinstance(workers, int), TypeError(
                f"'workers' expected a type 'int', got  {type(workers).__name__}."
            )
            assert workers >= 0, ValueError(
                f"'workers' must be positive or zero, got {workers}."
            )
        self.workers = workers

        self._compute()

    def is_unique(self, structure: Structure) -> bool:
        """Whether given structure is equivalent to any stored unique structure."""
        return any(self.matcher.fit(structure, struct) for struct in self.unique_structs)

    def _compute(self) -> None:
        # Remove eventual unphysical structures that could make the matcher throw an error
        unmatchables = list(filter(lambda struct: struct.volume < 1, self.structures))
        matchable_structs = list(filter(lambda struct: struct.volume >= 1, self.structures))

        # Group by stoichiometry
        formula_groups: list[list[Structure]] = group_compositions(
            matchable_structs, by="formula" # type: ignore
        )
        # Group by matching equivalence
        if self.workers is not None and self.workers == 0:
            groups: list[list[list[Structure]]] = [
                self.matcher.group_structures(
                    formula_grp
                ) for formula_grp in formula_groups
            ]
        else:
            groups: list[list[list[Structure]]] = process_map(
                self.matcher.group_structures,
                formula_groups,
                max_workers=self.workers,
                chunksize=min(10, len(formula_groups) // 100 + 1),
                desc="Matching structures for Unicity"
            )
        flattened_groups = list(itt.chain.from_iterable(groups))

        self._unique_structs = [group[0] for group in flattened_groups]
        self._duplicate_structs = list(itt.chain.from_iterable(
            [group[1:] for group in flattened_groups]
        ))
        self._unmatchable_structs = unmatchables

    @property
    def unique_structs(self) -> list[Structure]:
        """List of unique structures."""
        return self._unique_structs

    @property
    def duplicate_structs(self) -> list[Structure]:
        """List of duplicate structures."""
        return self._duplicate_structs

    @property
    def unmatchable_structs(self) -> list[Structure]:
        """List of structures that cannot be matched due to their unphysical volume."""
        return self._unmatchable_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Unique", self._unique_structs),
                ("Duplicate", self._duplicate_structs),
                ("Unmatchable", self._unmatchable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
        