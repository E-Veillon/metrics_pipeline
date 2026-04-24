"""Compute Novelty metric."""

import itertools as itt
from collections import OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from .metric_base import Metric
from .backend import group_compositions


class Novelty(Metric):
    """
    Compute Novelty metric.
    
    Definition
    ----------
    Structures are novel if they are not an equivalent representation of a reference structure.
    """
    def __init__(
        self,
        structures: list[Structure],
        ref_structs: list[Structure],
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5.0,
        workers: int | None = None,
        database_mode: bool = False
    ) -> None:
        """
        Compute Novelty metric.

        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate Novelty on.

        ref_structs: list[Structure]
            Known structures dataset to compare Novelty against.

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
            pass 0 to disable `process_map()` and do it sequentially. Unused if
            `database_mode` is set to True.

        database_mode: bool
            If set to True, automatic metric computation on given data is disabled,
            so you can provide an empty list to `structures` argument without error and a
            reference dataset of structures to `ref_structs` that will be stored to be compared
            against one structure at a time using the `is_novel()` method. Can be more efficient
            than automatic metric computation if the reference dataset is way larger than the
            number of structures to match. Defaults to False.
        """
        super().__init__(structures)
        self.ref_structs = ref_structs
        self.matcher = StructureMatcher(
            ltol=ltol, stol=stol, angle_tol=angle_tol, scale=False, attempt_supercell=True
        )

        if not database_mode:
            if workers is not None:
                if not isinstance(workers, int):
                    raise TypeError(
                        f"'workers' expected a type 'int', got  {type(workers).__name__}."
                )
                if workers < 0:
                    raise ValueError(
                        f"'workers' must be positive or zero, got {workers}."
                )
            self.workers = workers
            self._compute()

    def is_novel(self, structure: Structure) -> bool:
        """Whether given structure is equivalent to one of the initialized reference structures."""
        return any(self.matcher.fit(structure, ref_struct) for ref_struct in self.ref_structs)

    def _compute(self) -> None:
        for struct in self.structures:
            struct.properties["tmp_category"] = "computed"

        for struct in self.ref_structs:
            struct.properties["tmp_category"] = "reference"

        # Remove eventual unphysical structures that could make the matcher throw an error
        unmatchables = list(filter(lambda struct: struct.volume < 1, self.structures))
        matchable_structs = list(filter(lambda struct: struct.volume >= 1, self.structures))

        # Group by stoichiometry
        formula_groups: list[list[Structure]] = group_compositions(
            matchable_structs + self.ref_structs, by="formula" # type: ignore
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
                desc="Searching for novel structures"
            )
        flattened_groups = list(itt.chain.from_iterable(groups))
        self._novel_structs = list(itt.chain.from_iterable(
            [
                group for group in flattened_groups
                if all(struct.properties["tmp_category"] == "computed" for struct in group)
            ]
        ))
        self._known_structs = list(filter(
            lambda s: s.properties["tmp_category"] == "computed",
            itt.chain.from_iterable(
                [
                    group for group in flattened_groups
                    if any(struct.properties["tmp_category"] == "reference" for struct in group)
                ]
            )
        ))
        self._unmatchable_structs = unmatchables
        for struct in (self._novel_structs + self._known_structs + self._unmatchable_structs):
            struct.properties.pop("tmp_category")

    @property
    def novel_structs(self) -> list[Structure]:
        """List of novel structures."""
        return self._novel_structs

    @property
    def known_structs(self) -> list[Structure]:
        """List of structures equivalent to one in the reference dataset."""
        return self._known_structs

    @property
    def unmatchable_structs(self) -> list[Structure]:
        """List of structures that cannot be matched due to their unphysical volume."""
        return self._unmatchable_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Novel", self._novel_structs),
                ("Known", self._known_structs),
                ("Unmatchable", self._unmatchable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
