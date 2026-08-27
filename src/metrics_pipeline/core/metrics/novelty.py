"""Compute Novelty metric."""

import typing as tp
import itertools as itt
from collections import OrderedDict, defaultdict
import warnings

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from .metric_base import Metric, MetricsData
from metrics_pipeline.core.utils.genmat_data import GenMatStructure, StructureLike


class Novelty(Metric):
    """
    Compute Novelty metric.
    
    Definition
    ----------
    Structures are novel if they are not an equivalent representation of a reference structure.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
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
        structures: list[GenMatStructure]
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
            against one structure at a time using the `get_novelty()` or `is_novel()` methods.
            Can be more efficient than automatic metric computation if the reference dataset
            is way larger than the number of structures to match. Defaults to False.
        """
        super().__init__(structures, workers)
        self._matcher = StructureMatcher(
            ltol=ltol, stol=stol, angle_tol=angle_tol, scale=False, attempt_supercell=True
        )
        init_length = len(ref_structs)
        ref_structs = [struct for struct in ref_structs if struct.volume >= 1]
        if (num_unphysical_refs:=init_length - len(ref_structs)) > 0:
            warnings.warn(
                f"{type(self).__name__}: {num_unphysical_refs} reference structures having "
                "an unphysical volume < 1 angstrom were removed at initialization to avoid "
                "encountering bugs while matching."
            )
        self._ref_structs = self.sort_by_formula(ref_structs)

        if not database_mode:
            self._sorted_structs = self.sort_by_formula(structures)
            self._compute()

    @tp.overload
    def sort_by_formula(self, structures: list[Structure]) -> dict[str, list[Structure]]:
        ...
    @tp.overload
    def sort_by_formula(self, structures: list[GenMatStructure]) -> dict[str, list[GenMatStructure]]:
        ...
    def sort_by_formula(self, structures: list[StructureLike]) -> dict[str, list[StructureLike]]:
        """Sort structures by their reduced formula."""
        formulas = defaultdict(list)
        for structure in structures:
            formulas[structure.reduced_formula].append(structure)
        return formulas

    def is_novel(self, structure: StructureLike) -> bool:
        """
        Whether given structure is equivalent to one of the initialized reference structures.
        Compatible with standard `Structure` objects and return a simple bool instead of a
        complete `MetricsData` result object. Note that unphysical structures of volume < 1
        angstrom³ cannot be matched and are considered novel by default.
        """
        if structure.volume < 1 or structure.reduced_formula not in self._ref_structs:
            return True

        ref_structs = self._ref_structs[structure.reduced_formula]
        return any(
            self._matcher.fit(structure, ref_struct)
            for ref_struct in ref_structs
        )

    def get_novelty(self, structure: GenMatStructure) -> MetricsData:
        """
        Compute Novelty on a structure and get the results in a `MetricsData` object.
        """
        return MetricsData(
            structure,
            is_novel=self.is_novel(structure),
            is_unmatchable=structure.volume < 1,
            additional_data=self._get_metric_settings()
        )

    def _get_metric_settings(self) -> dict[str, tp.Any]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "ltol": self._matcher.ltol,
                "stol": self._matcher.stol,
                "angle_tol": self._matcher.angle_tol
            }
        }

    def _compute(self) -> None:
        self._computed_data: list[MetricsData] = []
        matched_structs: dict[str, list[GenMatStructure]] = defaultdict(list)

        # Only structures with formula in reference and with valid volumes are matched
        computed_extend = self._computed_data.extend
        computed_append = self._computed_data.append
        for formula, structs in self._sorted_structs.items():
            if formula not in self._ref_structs:
                computed_extend(
                    MetricsData(
                        struct,
                        is_novel=True,
                        additional_data=self._get_metric_settings()
                    ) for struct in structs
                )
            else:
                matched_append = matched_structs[formula].append
                for struct in structs:
                    if struct.volume < 1:
                        computed_append(
                            MetricsData(
                                struct,
                                is_unmatchable=True,
                                additional_data=self._get_metric_settings()
                            )
                        )
                    else:
                        matched_append(struct)

        # Group by stoichiometry
        formula_groups = (
            matched_structs[formula] + self._ref_structs[formula]
            for formula in matched_structs.keys()
        )
        # Group by matching equivalence
        num_groups = len(matched_structs)
        grouper = self._matcher.group_structures
        if self.workers == 0:
            groups: itt.chain[list[StructureLike]] = itt.chain.from_iterable(
                self._sequential_compute(
                    grouper, formula_groups, unit="systems matched",
                    lazy_data=True, n_elts=num_groups
                )
            )
        else:
            groups: itt.chain[list[StructureLike]] = itt.chain.from_iterable(
                self._parallel_compute(grouper, formula_groups)
            )

        for group in groups:
            if all(isinstance(struct, GenMatStructure) for struct in group):
                computed_extend(
                    MetricsData(
                        struct,
                        is_novel=True,
                        is_unmatchable=False,
                        additional_data=self._get_metric_settings()
                    ) for struct in group
                )
            else:
                computed_extend(
                    MetricsData(
                        struct,
                        is_novel=False,
                        is_unmatchable=False,
                        additional_data=self._get_metric_settings()
                    ) for struct in group if isinstance(struct, GenMatStructure)
                )

    @property
    def computed_data(self) -> list[MetricsData]:
        """List of all computed MetricsData objects."""
        return self._computed_data

    @property
    def novel_data(self) -> list[MetricsData]:
        """List of computed MetricsData containing a novel structure."""
        return [data for data in self.computed_data if data.is_novel]

    @property
    def known_data(self) -> list[MetricsData]:
        """
        List of computed MetricsData containing a structure equivalent
        to one in the reference dataset.
        """
        return [data for data in self.computed_data if data.is_novel is False]

    @property
    def unmatchable_data(self) -> list[MetricsData]:
        """
        List of MetricsData containing a structure that cannot be matched
        due to its unpysical volume.
        """
        return [data for data in self.computed_data if data.is_unmatchable]

    @property
    def novel_structs(self) -> list[GenMatStructure]:
        """List of novel structures."""
        return [data.typed_structure for data in self.computed_data if data.is_novel]

    @property
    def known_structs(self) -> list[GenMatStructure]:
        """List of structures equivalent to one in the reference dataset."""
        return [data.typed_structure for data in self.computed_data if data.is_novel is False]

    @property
    def unmatchable_structs(self) -> list[GenMatStructure]:
        """List of structures that cannot be matched due to their unphysical volume."""
        return [data.typed_structure for data in self.computed_data if data.is_unmatchable]

    @property
    def novel_names(self) -> list[str]:
        """List of names of novel structures."""
        return [data.typed_structure.name for data in self.computed_data if data.is_novel]

    @property
    def known_names(self) -> list[str]:
        """List of names of structures equivalent to one in the reference dataset."""
        return [data.typed_structure.name for data in self.computed_data if data.is_novel is False]

    @property
    def unmatchable_names(self) -> list[str]:
        """
        List of names of of structures that cannot be matched due to their unphysical volume.
        """
        return [data.typed_structure.name for data in self.computed_data if data.is_unmatchable]

    @property
    def novel_names_set(self) -> set[str]:
        """Set of names of novel structures."""
        return {data.typed_structure.name for data in self.computed_data if data.is_novel}

    @property
    def known_names_set(self) -> set[str]:
        """Set of names of structures equivalent to one in the reference dataset."""
        return {data.typed_structure.name for data in self.computed_data if data.is_novel is False}

    @property
    def unmatchable_names_set(self) -> set[str]:
        """Set of names of structures that cannot be matched due to their unphysical volume."""
        return {data.typed_structure.name for data in self.computed_data if data.is_unmatchable}

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[GenMatStructure]] = OrderedDict(
            [
                ("Novel", self.novel_structs),
                ("Known", self.known_structs),
                ("Unmatchable", self.unmatchable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
