"""Base class for implementing metrics classes."""

import typing as tp
from collections import OrderedDict
from abc import ABC, abstractmethod

import numpy as np

from pymatgen.core import Structure


class StructureDistribution(tp.Protocol):
    """
    Any callable taking a list of Structure objects and eventual keyword arguments
    and returning a numpy array representation of the structures distribution.
    """
    def __call__(self, structures: list[Structure], *args, **kwargs) -> np.ndarray:
        ...


class Metric(ABC):
    """Base class for implementing metrics classes. Do not call directly."""
    structures: list[Structure]

    def __init__(self, structures: list[Structure]) -> None:
        """Base class for implementing metrics classes."""
        assert isinstance(structures, list)
        for struct in structures:
            assert isinstance(struct, Structure)

        self.structures = structures

    def __len__(self) -> int:
        return len(self.structures)

    def __getitem__(self, idx: int) -> Structure:
        return self.structures[idx]

    @abstractmethod
    def _compute(self) -> None:
        """Compute the metric on stored structures and store result in dedicated attributes."""
        raise NotImplementedError

    @abstractmethod
    def write_result(self, filename: str) -> None:
        """Write a file with results of the metric computation."""
        raise NotImplementedError

    def _write_filter_metric_result(
        self, filename: str, subsets: OrderedDict[str, list[Structure]], verbose: bool = False
    ) -> None:
        """
        Write formatted results for filter metrics.

        Parameters
        ----------
        filename: str
            Path to write the result.
 
        subsets: OrderedDict[str, list[Structure]]
            Ordered dictionary of each subset of structures to show in the report.
            The keys will be used to name each category. Should contain at least
            2 subsets of structures with respect to whether they passed or failed the filter.

        verbose: bool
            Whether to list the structures in each subset. Defaults to False.
        """
        text_lines = [f"===== {self.__class__.__name__} Results ====="]
        text_lines.append(f"Total structures:  {len(self)}")

        for name, subset in subsets.items():
            text_lines.append(f"{name} structures: {len(subset)}")

        if verbose:
            for name, subset in subsets.items():
                text_lines.append(f"\nList of {name} structures:")
                text_lines.extend(list(map(str, subset)))

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text_lines))
