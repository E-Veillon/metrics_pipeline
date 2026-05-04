"""Compute Unicity metric."""

import typing as tp
import itertools as itt
from collections import OrderedDict, defaultdict
from collections.abc import Callable

from pymatgen.analysis.structure_matcher import StructureMatcher

from .metric_base import Metric, MetricsData
from src.utils.genmat_data import GenMatStructure, StructureLike


class Unicity(Metric):
    """
    Compute Unicity metric.
    
    Definition
    ----------
    Structures are unique if they are not an equivalent representation of
    a previous structure in the same dataset.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5.0,
        workers: int | None = None
    ) -> None:
        """
        Compute Unicity metric.

        Parameters
        ----------
        structures: list[GenMatStructure]
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
        super().__init__(structures, workers)
        self._matcher = StructureMatcher(
            ltol=ltol, stol=stol, angle_tol=angle_tol, scale=False, attempt_supercell=True
        )
        self._compute()

    def is_unique(self, structure: StructureLike) -> bool:
        """
        Compute unicity of a structure compared to stored unique structures.
        Compatible with standard `Structure` objects, but returns a bool instead
        of a computed `MetricsData` object.

        Parameters
        ----------
        structure: Structure
            The structure to compute Unicity on.

        Returns
        -------
        bool
            Whether the structure is unique compared to stored dataset.
        """
        # Check matchability
        if structure.volume < 1:
            return True
        # Build the list of unique and matchable structures with same reduced formula
        formula = structure.reduced_formula
        predicate: Callable[[MetricsData], bool] = lambda d: (
            not d.is_unmatchable and formula == d.structure.reduced_formula
        )
        ref_structs = [data.typed_structure for data in self.unique_data if predicate(data)]
        # If no structure with same reduced formula, it's unique
        if not ref_structs:
            return True
        # Match with structures with same reduced formula
        return not any(self._matcher.fit(structure, struct) for struct in ref_structs)

    def get_unicity(self, structure: GenMatStructure, add_to_data: bool = False) -> MetricsData:
        """
        Compute unicity of a structure compared to stored unique structures.
        Return a computed `MetricsData` containing the structure and metric results.

        Parameters
        ----------
        structure: GenMatStructure
            The structure to compute Unicity on.

        add_to_data: bool
            Whether to add the newly computed `MetricsData` to the metric dataset.
            Defaults to False.

        Returns
        -------
        MetricsData
            Computed metric results.
        """
        data = MetricsData(
            structure,
            is_unique=self.is_unique(structure),
            is_unmatchable=structure.volume < 1,
            additional_data=self._get_metric_settings()
        )
        if add_to_data:
            self._computed_data.append(data)
        return data

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

        for struct in self.structures:
            if struct.volume < 1:
                # Structures with unphysical volume are removed from computation
                self._computed_data.append(
                    MetricsData(
                        struct,
                        is_unmatchable=True,
                        additional_data=self._get_metric_settings()
                    )
                )
            else:
                # Sort matchable structures by formula
                matched_structs[struct.reduced_formula].append(struct)

        # Group by matching equivalence
        if self.workers == 0:
            groups: itt.chain[list[GenMatStructure]] = itt.chain.from_iterable(
                self._sequential_compute(
                    self._matcher.group_structures, matched_structs.values()
                )
            )
        else:
            groups: itt.chain[list[GenMatStructure]] = itt.chain.from_iterable(
                self._parallel_compute(
                    self._matcher.group_structures, matched_structs.values()
                )
            )
        for group in groups:
            self._computed_data.append(
                MetricsData(
                    group[0],
                    is_unique=True,
                    is_unmatchable=False,
                    additional_data=self._get_metric_settings()
                )
            )
            if len(group) > 1:
                self._computed_data.extend(
                    MetricsData(
                        struct,
                        is_unique=False,
                        is_unmatchable=False,
                        additional_data=self._get_metric_settings()
                    ) for struct in group[1:]
                )

    @property
    def computed_data(self) -> list[MetricsData]:
        """List of all computed MetricsData objects."""
        return self._computed_data

    @property
    def unique_data(self) -> list[MetricsData]:
        """List of computed MetricsData containing a unique structure."""
        return [data for data in self.computed_data if data.is_unique]

    @property
    def duplicate_data(self) -> list[MetricsData]:
        """List of computed MetricsData containing duplicate structures."""
        return [data for data in self.computed_data if data.is_unique is False]

    @property
    def unmatchable_data(self) -> list[MetricsData]:
        """
        List of computed MetricsData containing a structure that cannot be matched
        due to its unphysical volume.
        """
        return [data for data in self.computed_data if data.is_unmatchable]

    @property
    def unique_structs(self) -> list[GenMatStructure]:
        """List of unique structures."""
        return [data.typed_structure for data in self.computed_data if data.is_unique]

    @property
    def duplicate_structs(self) -> list[GenMatStructure]:
        """List of duplicate structures."""
        return [data.typed_structure for data in self.computed_data if data.is_unique is False]

    @property
    def unmatchable_structs(self) -> list[GenMatStructure]:
        """List of structures that cannot be matched due to their unphysical volume."""
        return [data.typed_structure for data in self.computed_data if data.is_unmatchable]

    @property
    def unique_names(self) -> list[str]:
        """List of names of unique structures."""
        return [data.typed_structure.name for data in self.computed_data if data.is_unique]

    @property
    def duplicate_names(self) -> list[str]:
        """List of names of duplicate structures."""
        return [data.typed_structure.name for data in self.computed_data if data.is_unique is False]

    @property
    def unmatchable_names(self) -> list[str]:
        """List of names of structures that cannot be matched due to their unphysical volume."""
        return [data.typed_structure.name for data in self.computed_data if data.is_unmatchable]

    @property
    def unique_names_set(self) -> set[str]:
        """Set of names of unique structures."""
        return {data.typed_structure.name for data in self.computed_data if data.is_unique}

    @property
    def duplicate_names_set(self) -> set[str]:
        """Set of names of duplicate structures."""
        return {data.typed_structure.name for data in self.computed_data if data.is_unique is False}

    @property
    def unmatchable_names_set(self) -> set[str]:
        """Set of names of structures that cannot be matched due to their unphysical volume."""
        return {data.typed_structure.name for data in self.computed_data if data.is_unmatchable}

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[GenMatStructure]] = OrderedDict(
            [
                ("Unique", self.unique_structs),
                ("Duplicate", self.duplicate_structs),
                ("Unmatchable", self.unmatchable_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
        