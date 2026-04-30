"""Base class for implementing metrics classes."""

import typing as tp
import typing_extensions as tpe
from collections import OrderedDict
from collections.abc import Callable
from abc import ABC, abstractmethod
from dataclasses import dataclass, asdict, field

import numpy as np
from tqdm.contrib.concurrent import process_map

from src.utils import VisualIterator, check_type, check_num_value
from src.utils.genmat_data import GenMatPDEntry, GenMatStructure, StructureLike


class StructureDistribution(tp.Protocol):
    """
    Any callable taking a list of Structure objects and eventual keyword arguments
    and returning a numpy array representation of the structures distribution.
    """
    def __call__(self, structures: list[StructureLike], *args, **kwargs) -> np.ndarray:
        ...


class Metric(ABC):
    """Base class for implementing metrics classes. Do not call directly."""
    __metric_properties__: frozenset[str] = frozenset()
    structures: list[GenMatStructure]

    def __init_subclass__(cls, *args, **kwargs) -> None:
        super().__init_subclass__(*args, **kwargs)
        # Store all property method names in a custom magic class attribute
        properties = {
            name for name, value in cls.__dict__.items() if isinstance(value, property)
        }
        # If the Metric class inherits from another one, also store base class properties
        for base in cls.mro():
            properties.update(
                {name for name, value in base.__dict__.items() if isinstance(value, property)}
            )
        cls.__metric_properties__ = frozenset(properties)

    def __init__(self, structures: list[GenMatStructure], workers: int | None = None) -> None:
        """Base class for implementing metrics classes."""
        assert isinstance(structures, list)
        for idx, struct in enumerate(structures):
            check_type(struct, f"structures[{idx}]", (GenMatStructure,))
        check_type(workers, "workers", (int, type(None)))
        if workers is not None:
            check_num_value(workers, "workers", ">=", 0)

        self.structures = structures
        self.workers = workers
        self.chunksize = min(len(structures) // 100 + 1, 10)

    def __len__(self) -> int:
        return len(self.structures)

    def __getitem__(self, idx: int) -> GenMatStructure:
        return self.structures[idx]

    def __setattr__(self, name: str, value: tp.Any) -> None:
        raise AttributeError(f"Cannot assign attributes ({type(self).__name__!r} is immutable).")

    def __delattr__(self, name: str) -> None:
        raise AttributeError(f"Cannot delete attributes ({type(self).__name__!r} is immutable).")

    def _sequential_compute(
        self,
        func: Callable,
        data: tp.Iterable[tp.Any],
        unpack: bool = False,
        desc: str | None = None,
        unit: str | None = None,
        lazy_data: bool = False,
        n_elts: int | None = None
    ) -> tp.Generator[tp.Any, None, None]:
        """
        Initialize a `VisualIterator` and compute the results lazily.
        
        Parameters
        ----------
        func: Callable
            Function to compute on the data.

        data: Iterable[Any]
            The source data to compute on.

        unpack: bool
            Whether to try to unpack `data` items when passing them to `func`.
            Defaults to False.

        desc: str
            Task description to print on the console. Defaults to 'Computing a metric'.

        unit: str
            Suffix to print after the progression indicator. Defaults to 'computed'.

        lazy_data: bool
            Whether `data` is a big iterator that should not be expanded into memory.
            Defaults to False.

        n_elts: int, optional
            Length of `data`. Only used if `lazy_data` is set to `True` to print a total
            length and percentage on progression indicator.
        """
        if desc is None:
            desc = f"Computing {type(self).__name__} metric"

        if unit is None:
            unit = "computed"

        if lazy_data:
            iterator = VisualIterator.from_big_iterator(
                iter(data), n_elts, desc=desc, unit=unit, percent=True
            )
        else:
            iterator = VisualIterator(data, desc, unit, percent=True)

        for item in iterator:
            if unpack:
                yield func(*item)
            else:
                yield func(item)

    def _parallel_compute(
        self,
        func: Callable,
        data: tp.Iterable[tp.Any],
        unpack: bool = False,
        desc: str | None = None
    ) -> list[tp.Any]:
        """
        Compute using parallel processes with initialized `workers` value.

        Parameters
        ----------
        func: Callable
            Function to compute on the data.

        data: Iterable[Any]
            The source data to compute on.

        unpack: bool
            Whether to try to unpack `data` items when passing them to `func`.
            Defaults to False.

        desc: str
            Task description to print on the console. Defaults to 'Computing a metric'.
        """
        if desc is None:
            desc = f"Computing {type(self).__name__} metric"

        if unpack:
            return process_map(
            func, *data, max_workers=self.workers, chunksize=self.chunksize, desc=desc
        )
        return process_map(
            func, data, max_workers=self.workers, chunksize=self.chunksize, desc=desc
        )

    @abstractmethod
    def _get_metric_settings(self) -> dict[str, tp.Any]:
        """Get a dict of initialized parameters for this metric."""
        raise NotImplementedError

    @abstractmethod
    def _compute(self) -> None:
        """Compute the metric on stored structures and store result in dedicated attributes."""
        raise NotImplementedError

    @abstractmethod
    def write_result(self, filename: str) -> None:
        """Write a file with results of the metric computation."""
        raise NotImplementedError

    def _write_filter_metric_result(
        self, filename: str, subsets: OrderedDict[str, list[GenMatStructure]], verbose: bool = False
    ) -> None:
        """
        Write formatted results for filter metrics.

        Parameters
        ----------
        filename: str
            Path to write the result.
 
        subsets: OrderedDict[str, list[GenMatStructure]]
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
                text_lines.extend(f"- {struct.name}" for struct in subset)

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(text_lines))


class MetricNotComputedError(AttributeError):
    """Metric was not computed, the property is not defined."""
    ...


class NoneMetric:
    """
    Placeholder metric class to substitute non-computed metrics.
    Raises `MetricNotComputedError` whenever a wrapped property is accessed.
    """
    def __init__(self, wrapped_metric: type[Metric]) -> None:
        """
        Placeholder metric class to substitute non-computed metrics.
        Raises `MetricNotComputedError` whenever a wrapped property is accessed.
        """
        check_type(wrapped_metric, "wrapped_metric", (type,))
        if not issubclass(wrapped_metric, Metric):
            raise TypeError(
                f"{type(self).__name__}: wrapped class must inherit from "
                f"the 'Metric' Abstract Base Class, got {wrapped_metric.__name__!r}."
            )
        self._wrapped_metric_name = wrapped_metric.__name__
        self.__metric_properties__: frozenset = frozenset(wrapped_metric.__metric_properties__)

    def __getattr__(self, name: str) -> tp.Any:
        if name in self.__metric_properties__:
            raise MetricNotComputedError(
                f"Metric {self._wrapped_metric_name!r} was not computed, "
                f"property {name!r} is not available."
            )
        raise AttributeError(f"{self._wrapped_metric_name!r} has no attribute {name!r}.")

    def __setattr__(self, name: str, value: tp.Any) -> None:
        raise AttributeError(f"Cannot assign attributes ({type(self).__name__!r} is immutable).")
    
    def __delattr__(self, name: str) -> None:
        raise AttributeError(f"Cannot delete attributes ({type(self).__name__!r} is immutable).")

    def __bool__(self) -> bool:
        return False

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self._wrapped_metric_name})"

    def __dir__(self) -> list[str]:
        return sorted(self.__metric_properties__)


@dataclass
class MetricsData:
    """
    Dataclass to store all metrics data associated with a tested structure.

    Attributes
    ----------
    structure: GenMatStructure
        The underlying GenMatStructure object.
    is_valid: bool | None
        Whether the structure is valid, i.e. no pair of atoms are less than 0.5 Angstrom apart.
    is_viable: bool | None
        Whether the structure is viable, i.e. no pair of atoms are closer than the sum of their
        atomic radii within a tolerance and no lattice parameter is unphysical.
    is_symmetric: bool | None
        Whether the structure is symmetric, i.e. of higher symmetry than triclinic system
        (spacegroups 1 and 2).
    is_element_metastable: bool | None
        Whether the structure has negative formation energy relative to reference elements.
    is_metastable: bool | None
        Whether the structure energy above hull is below the metastability threshold.
    is_stable: bool | None
        Whether the structure energy above hull is negative or zero.
    is_unique: bool | None
        Whether the structure is unique, i.e. not a duplicate of another structure
        in tested set.
    is_novel: bool | None
        Whether the structure is novel, i.e. not a duplicate of a structure
        in a reference dataset.
    is_unmatchable: bool | None
        Whether the structure is unmatchable, i.e. has an unphysically small volume
        that makes it impossible to use StructureMatcher on it for Unicity or Novelty.
    rmsd: float | None
        The RMSD of the structure with respect to its relaxed state, in Angstrom.
    additional_data: dict
        A dictionary to store any additional data related to metrics computations.
    """
    METRIC_SLOTS = (
        "is_valid", "is_viable", "is_symmetric",
        "is_element_metastable", "is_metastable", "is_stable",
        "is_unique", "is_novel", "is_unmatchable",
        "rmsd"
    )
    __slots__ = ("structure", *METRIC_SLOTS, "additional_data")

    structure: GenMatStructure
    is_valid: bool | None = None
    is_viable: bool | None = None
    is_symmetric: bool | None = None
    is_element_metastable: bool | None = None
    is_metastable: bool | None = None
    is_stable: bool | None = None
    is_unique: bool | None = None
    is_novel: bool | None = None
    is_unmatchable: bool | None = None
    rmsd: float | None = None
    additional_data: dict = field(default_factory=dict)

    def __post_init__(self):
        check_type(self.structure, "structure", (GenMatStructure,))
        check_type(self.rmsd, "rmsd", (float, type(None)))
        if self.rmsd is not None:
            check_num_value(self.rmsd, "rmsd", ">=", 0)
        check_type(self.additional_data, "additional_data", (dict,))
        if self.is_unmatchable:
            self.is_unique = True
            self.is_novel = True

    def update(self, other: tpe.Self) -> None:
        """Update the current MetricsData attributes using another MetricsData object."""
        check_type(other, "other", (type(self),))
        if self.name != other.name:
            raise ValueError(
                f"Cannot update MetricsData with different names: "
                f"{self.name!r} and {other.name!r}."
            )
        for field_name in self.METRIC_SLOTS:
            if (new_value := getattr(other, field_name, None)) is not None:
                setattr(self, field_name, new_value)

    def copy(self) -> tpe.Self:
        """Create a new instance with same structure and metrics values."""
        cls = type(self)
        new = cls(self.structure)

        for metric_slot in self.METRIC_SLOTS:
            setattr(new, metric_slot, getattr(self, metric_slot))
        
        new.additional_data |= self.additional_data

        return new

    def __or__(self, other: tpe.Self) -> tpe.Self:
        if not isinstance(other, type(self)):
            return NotImplemented
        
        new = self.copy()
        new.update(other)
        return new

    def __ior__(self, other: tpe.Self)-> tpe.Self:
        if not isinstance(other, type(self)):
            return NotImplemented
        
        self.update(other)
        return self

    @property
    def typed_structure(self) -> GenMatStructure:
        """
        The internal structure, equivalent to `structure` attribute.
        Only redefined to get more explicit type checking.
        """
        return self.structure

    @property
    def name(self) -> str:
        """The name associated with the structure."""
        return self.typed_structure.name

    @property
    def entry(self) -> GenMatPDEntry:
        """The `GenMatPDEntry` object associated with the structure."""
        return self.typed_structure.entry

    @property
    def energy_per_atom(self) -> float | None:
        """The energy per atom of the structure, in eV/atom."""
        return self.typed_structure.energy_per_atom

    def as_dict(self) -> dict[str, tp.Any]:
        """Get a dictionary representation of the MetricsData object."""
        dct = asdict(self)
        dct["structure"] = self.typed_structure.as_dict()
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> tpe.Self:
        """Create a MetricsData object from a dictionary representation."""
        return cls(
            structure=GenMatStructure.from_dict(dct["structure"]),
            is_valid=dct.get("is_valid"),
            is_viable=dct.get("is_viable"),
            is_symmetric=dct.get("is_symmetric"),
            is_element_metastable=dct.get("is_element_metastable"),
            is_metastable=dct.get("is_metastable"),
            is_stable=dct.get("is_stable"),
            is_unique=dct.get("is_unique"),
            is_novel=dct.get("is_novel"),
            is_unmatchable=dct.get("is_unmatchable"),
            rmsd=dct.get("rmsd"),
            additional_data=dct.get("additional_data", {})
        )
