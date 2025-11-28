"""Base class for implementing metrics classes."""

from abc import ABC, abstractmethod
from pymatgen.core import Structure


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
    def write_result(self, filename: str, verbose: bool = False) -> None:
        """Write a file with results of the metric computation."""
        raise NotImplementedError

    def _write_filter_metric_result(
        self,
        attr_ok: list[Structure], attr_not_ok: list[Structure],
        filename: str, verbose: bool = False
    ) -> None:
        """
        Write formatted results for filter metrics.
        
        Parameters
        ----------
        attr_ok: list[Structure]
            Attribute containing structures that passed the filter.

        attr_not_ok: list[Structure]
            Attribute containing structures that did not pass the filter.

        filename: str
            Path to write the result.

        verbose: bool
            Whether to write the lists of all structures that passed or failed the filter.
            Defaults to False.
        """
        text = f"===== {self.__class__.__name__} Result ====="
        text += f"Total structures:  {len(self.structures)}"
        text += f"Passed structures: {len(attr_ok)}"
        text += f"Failed structures: {len(attr_not_ok)}"

        if verbose:
            passed_list = "\n".join(list(map(str, attr_ok)))
            text += f"List of structures that passed:\n{passed_list}\n"
            failed_list = "\n".join(list(map(str, attr_not_ok)))
            text += f"List of structures that failed:\n{failed_list}"

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text))

    def _write_similarity_metric_result(self) -> None:
        """
        Write formatted results for similarity metrics.
        """
        raise NotImplementedError
