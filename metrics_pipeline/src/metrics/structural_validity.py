"""Compute Structural Validity metric."""

from collections import OrderedDict
from typing import Any

from .metric_base import Metric, MetricsData
from src.utils import GenMatStructure, StructureLike


class StructValidity(Metric):
    """
    Compute Structural Validity metric.
    
    Definition
    ----------
    A structure is valid if none of its atoms are closer than a threshold, generally 0.5 angstroms.

    Reference
    ---------
    Court, C. J., Yildirim, B., Jain, A., & Cole, J. M. (2020).
    3-D inorganic crystal structure generation and property prediction via representation learning.
    Journal of Chemical Information and Modeling, 60(10), 4518-4535.
    """
    def __init__(
        self,
        structures: list[GenMatStructure],
        valid_tol: float = 0.5,
        workers: int | None = None
    ) -> None:
        """
        Compute Structural Validity metric.

        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Structural Validity on.

        valid_tol: float
           Structural Validity threshold for atoms distance.
           Defaults to 0.5 angstroms.

        workers: int, optional
            Number of parallel processes to spawn for high throughput structure matching.
            If not given, will use default of `tqdm.contrib.concurrent.process_map()`.
            pass 0 to disable `process_map()` and do it sequentially.
        """
        super().__init__(structures, workers)
        self._valid_tol = valid_tol

        self._compute()

    def is_valid(self, structure: StructureLike) -> bool:
        """
        Compute Validity of a structure relative to initialized tolerance.
        Compatible with standard `Structure` objects, but only return a bool
        instead of a computed `MetricsData` object with metric results.
        """
        return structure.is_valid(tol=self._valid_tol)

    def get_validity(self, structure: GenMatStructure) -> MetricsData:
        """
        Compute Validity of a structure relative to initialized tolerance,
        and return a computed MetricsData object of the metric results.
        """
        return MetricsData(
            structure,
            is_valid=self.is_valid(structure),
            additional_data=self._get_metric_settings()
        )

    def _get_metric_settings(self) -> dict[str, Any]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "valid_tol": self._valid_tol
            }
        }

    def _compute(self) -> None:
        if self.workers == 0:
            self._computed_data = list(self._sequential_compute(
                    self.get_validity, self.structures
                ))
        else:
            self._computed_data = self._parallel_compute(
                self.get_validity, self.structures
            )

    @property
    def computed_data(self) -> list[MetricsData]:
        """List of all computed MetricsData objects."""
        return self._computed_data

    @property
    def valid_data(self) -> list[MetricsData]:
        """List of computed MetricsData objects containing valid structures."""
        return [data for data in self.computed_data if data.is_valid]

    @property
    def invalid_data(self) -> list[MetricsData]:
        """List of computed MetricsData objects containing invalid structures."""
        return [data for data in self.computed_data if data.is_valid is False]

    @property
    def valid_structs(self) -> list[GenMatStructure]:
        """List of valid structures."""
        return [data.typed_structure for data in self.computed_data if data.is_valid]

    @property
    def invalid_structs(self) -> list[GenMatStructure]:
        """List of invalid structures."""
        return [data.typed_structure for data in self.computed_data if data.is_valid is False]

    @property
    def valid_names(self) -> list[str]:
        """List of names of valid structures."""
        return [data.typed_structure.name for data in self.computed_data if data.is_valid]
    
    @property
    def invalid_names(self) -> list[str]:
        """List of names of invalid structures."""
        return [data.typed_structure.name for data in self.computed_data if data.is_valid is False]

    @property
    def valid_names_set(self) -> set[str]:
        """Set of names of valid structures."""
        return {data.typed_structure.name for data in self.computed_data if data.is_valid}
    
    @property
    def invalid_names_set(self) -> set[str]:
        """Set of names of invalid structures."""
        return {data.typed_structure.name for data in self.computed_data if data.is_valid is False}

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[GenMatStructure]] = OrderedDict(
            [
                ("Valid", self.valid_structs),
                ("Invalid", self.invalid_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)