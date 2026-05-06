"""Compute the Stability, Unicity, Novelty compound metric (S.U.N.)."""

import typing as tp
from collections import OrderedDict

from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PDEntry

from .metric_base import Metric, NoneMetric, MetricsData
from .stability import Stability
from .unicity import Unicity
from .novelty import Novelty
from src.utils import GenMatStructure


# TODO: Tester le comportement des propriétés internes avant de lancer en production!
class SUN(Metric):
    """
    Compute the Stability, Unicity, Novelty compound metric (S.U.N.).
    
    Definition
    ----------
    S.U.N. structures are high quality structures passing the three metrics of
    Stability, Unicity and Novelty (see their respective modules for more details on
    the definition of each of these metrics).

    Reference
    ---------
    Claudio Zeni et al. “A generative model for inorganic materials design”.
    In: Nature 639.8055 (2025), pp. 624–632.
    """
    _ref_structs: list[Structure]
    _compute_stability: bool
    _compute_unicity: bool
    _compute_novelty: bool
    _stable_tol: float
    _ltol: float
    _stol: float
    _angle_tol: float
    _stability: Stability | NoneMetric
    _unicity: Unicity | NoneMetric
    _novelty: Novelty | NoneMetric

    def __init__(
        self,
        structures: list[GenMatStructure],
        ref_structs: list[Structure],
        compute_stability: bool = True,
        compute_unicity: bool = True,
        compute_novelty: bool = True,
        workers: int | None = None,
        **kwargs
    ) -> None:
        """
        Compute the Stability, Unicity, Novelty compound metric (S.U.N.).
        
        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate S.U.N. metrics on.

        ref_structs: list[Structure]
            Structures of known formation energy to use as references.

        compute_stability: bool
            Whether to compute the Stability metric as part of S.U.N. computation.
            Defaults to True.

        compute_unicity: bool
            Whether to compute the Unicity metric as part of S.U.N. computation.
            Defaults to True.

        compute_novelty: bool
            Whether to compute the Novelty metric as part of S.U.N. computation.
            Defaults to True.

        workers: int, optional
            Number of parallel processes to spawn for high throughput structure matching.
            If not given, will use default of `tqdm.contrib.concurrent.process_map()`.
            pass 0 to disable `process_map()` and do it sequentially.

        kwargs: Any
            Additionnal arguments to pass to the respective metrics.
        
        Notes
        -----
        - Given structures (both computed and references) must have their total energy
        in eV stored in their properties under the 'energy' key for Stability computation.

        - Structures having an unphysical volume < 1 angstrom³ can cause a bug if passed
        to the pymatgen's StructureMatcher used for Unicity and Novelty. Hence, they are
        checked and removed from those computations beforehand and labeled as 'unmatchable'.
        Due to their very unlikely shape, unmatchable structures are all assumed to be unique
        and novel by default, but are not suitable for any kind of chemistry.

        - If a metric is deactivated, all structures are assumed to pass it by default.
        If neither Unicity nor Novelty is computed, all structures are assumed to be matchable.
        """
        super().__init__(structures, workers)
        self._ref_structs = ref_structs
        self._compute_stability = compute_stability
        self._compute_unicity = compute_unicity
        self._compute_novelty = compute_novelty
        self._stable_tol = kwargs.pop("stable_tol", 0.1)
        self._ltol = kwargs.pop("ltol", 0.2)
        self._stol = kwargs.pop("stol", 0.3)
        self._angle_tol = kwargs.pop("angle_tol", 5.0)
        self._database_mode = kwargs.pop("database_mode", False)

        self._compute()

    def _get_ref_entries(self) -> list[PDEntry]:
        """Convert reference structures to `PDEntry` objects."""
        return [
            PDEntry(s.composition, s.properties["energy"], s.properties.get("header"))
            for s in self._ref_structs
        ]

    def _get_metric_settings(self) -> dict[str, tp.Any]:
        """Get a dict of initialized parameters for this metric."""
        stability_settings = self._stability._get_metric_settings() if self._compute_stability else {}
        unicity_settings = self._unicity._get_metric_settings() if self._compute_unicity else {}
        novelty_settings = self._novelty._get_metric_settings() if self._compute_novelty else {}
        return stability_settings | unicity_settings | novelty_settings

    def _compute(self) -> None:
        # Remove unmatchable structures directly to avoid internal bugs
        if self._compute_unicity or self._compute_novelty:
            matchable_structs = []
            unmatchable_names = set()
            for struct in self.structures:
                if struct.volume < 1:
                    unmatchable_names.add(struct.name)
                else:
                    matchable_structs.append(struct)
        else:
            matchable_structs = self.structures
            unmatchable_names = set()

        self._stability = Stability(
            self.structures, self._get_ref_entries(), self._stable_tol, workers=self.workers
        ) if self._compute_stability else NoneMetric(Stability)

        metastable_names = self._stability.metastable_names_set if self._compute_stability else set()
        stable_names = self._stability.stable_names_set if self._compute_stability else set()

        self._unicity = Unicity(
            matchable_structs, self._ltol, self._stol, self._angle_tol, self.workers
        ) if self._compute_unicity else NoneMetric(Unicity)

        unique_names = self._unicity.unique_names_set if self._compute_unicity else set()

        if self._compute_novelty and self._database_mode:
            self._novelty = Novelty(
                [], self._ref_structs,
                self._ltol, self._stol, self._angle_tol, self.workers, self._database_mode
            )
            novel_names = {
                struct.name for struct in matchable_structs if self._novelty.is_novel(struct)
            }
        else:
            self._novelty = Novelty(
                matchable_structs, self._ref_structs,
                self._ltol, self._stol, self._angle_tol, self.workers, self._database_mode
            ) if self._compute_novelty else NoneMetric(Novelty)
            novel_names = self._novelty.novel_names_set if self._compute_novelty else set()

        self._computed_data: list[MetricsData] = [
            MetricsData(
                struct,
                is_metastable=(
                    struct.name in metastable_names
                    if self._compute_stability else None
                ),
                is_stable=(
                    struct.name in stable_names
                    if self._compute_stability else None
                ),
                is_unique=(
                    struct.name in unique_names
                    if self._compute_unicity else None
                ),
                is_novel=(
                    struct.name in novel_names
                    if self._compute_novelty else None
                ),
                is_unmatchable=(
                    struct.name in unmatchable_names
                    if self._compute_unicity or self._compute_novelty else None
                ),
                additional_data=self._get_metric_settings()
            ) for struct in self.structures
        ]

    @property
    def computed_data(self) -> list[MetricsData]:
        """List of all computed `MetricsData` objects."""
        return self._computed_data

    @property
    def unmatchable_data(self) -> list[MetricsData]:
        """
        List of computed `MetricsData` containing a structure that cannot be matched
        due to its unphysical volume.
        
        Notes
        -----
        When computing Unicity and/or Novelty through this class, unmatchable structures
        are removed before calling internal metrics to avoid conflicts and double
        computation, thus the `unmatchable_data` property from internal metrics is always empty,
        and this global property is the only populated one.
        """
        if self._compute_unicity or self._compute_novelty:
            return [data for data in self.computed_data if data.is_unmatchable]
        return NoneMetric(type(self)).unmatchable_data

    @property
    def unmatchable_structs(self) -> list[GenMatStructure]:
        """
        List of structures that cannot be matched due to their unphysical volume.
        
        Notes
        -----
        When computing Unicity and/or Novelty through this class, unmatchable structures
        are removed before calling internal metrics to avoid conflicts and double
        computation, thus the `unmatchable_structs` property from internal metrics is always empty,
        and this global property is the only populated one.
        """
        if self._compute_unicity or self._compute_novelty:
            return [data.typed_structure for data in self.computed_data if data.is_unmatchable]
        return NoneMetric(type(self)).unmatchable_structs

    @property
    def unmatchable_names(self) -> list[str]:
        """
        List of names of structures that cannot be matched due to their unphysical volume.
        
        Notes
        -----
        When computing Unicity and/or Novelty through this class, unmatchable structures
        are removed before calling internal metrics to avoid conflicts and double
        computation, thus the `unmatchable_names` property from internal metrics is always empty,
        and this global property is the only populated one.
        """
        if self._compute_unicity or self._compute_novelty:
            return [data.typed_structure.name for data in self.computed_data if data.is_unmatchable]
        return NoneMetric(type(self)).unmatchable_names

    @property
    def unmatchable_names_set(self) -> set[str]:
        """
        Set of names of structures that cannot be matched due to their unphysical volume.
        
        Notes
        -----
        When computing Unicity and/or Novelty through this class, unmatchable structures
        are removed before calling internal metrics to avoid conflicts and double
        computation, thus the `unmatchable_names_set` property from internal metrics is always empty,
        and this global property is the only populated one.
        """
        if self._compute_unicity or self._compute_novelty:
            return {data.typed_structure.name for data in self.computed_data if data.is_unmatchable}
        return NoneMetric(type(self)).unmatchable_names_set

    @property
    def stability(self) -> Stability | NoneMetric:
        """
        The computed Stability metric object, or `NoneMetric` placeholder if not computed.
        """
        return self._stability

    @property
    def unicity(self) -> Unicity | NoneMetric:
        """
        The computed Unicity metric object, or `NoneMetric` placeholder if not computed.
        """
        return self._unicity

    @property
    def novelty(self) -> Novelty | NoneMetric:
        """
        The computed Novelty metric object, or `NoneMetric` placeholder if not computed.
        """
        return self._novelty

    def get_computed_subset(
        self,
        stable: bool | None = None,
        unique: bool | None = None,
        novel: bool | None = None,
        unmatchable: bool | None = None
    ) -> list[MetricsData]:
        """
        Get a subset of computed data according to given metrics flags.

        Parameters
        ----------
        stable: bool or None
            - If set to True, only return data for structures that passed the Stability metric.
            - If set to False, only return data for structures that failed the Stability metric.
            - If set to None, return data for structures indifferently for this metric.
            Defaults to None.

        unique: bool or None
            - If set to True, only return data for structures that passed the Unicity metric.
            - If set to False, only return data for structures that failed the Unicity metric.
            - If set to None, return data for structures indifferently for this metric.
            Defaults to None.

        novel: bool or None
            - If set to True, only return data for structures that passed the Novelty metric.
            - If set to False, only return data for structures that failed the Novelty metric.
            - If set to None, return data for structures indifferently for this metric.
            Defaults to None.

        unmatchable: bool or None
            - If set to True, only return data for structures that could not be matched.
            - If set to False, only return data for structures that could be matched.
            - If set to None, return data for structures indifferently for this metric.
            Defaults to None.

        Returns
        -------
        list[MetricsData]
            The list of computed `MetricsData` with the queried metrics properties.

        Notes
        -----
        - Passing `None` to all metrics (default) is equivalent to calling the
        `computed_data` property.
        """
        if stable is None and unique is None and novel is None and unmatchable is None:
            return self.computed_data

        return list(
            data for data in filter(
                lambda s: (
                    (True if stable is None else s.is_stable is stable) and
                    (True if unique is None else s.is_unique is unique) and
                    (True if novel is None else s.is_novel is novel) and
                    (True if unmatchable is None else s.is_unmatchable is unmatchable)
                ), self.computed_data
            )
        )
    
    def get_structs_subset(
        self,
        stable: bool | None = None,
        unique: bool | None = None,
        novel: bool | None = None,
        unmatchable: bool | None = None
    ) -> list[GenMatStructure]:
        """
        Get the subset of structures according to given metrics flags.
        Behaves just like the `get_computed_subset()` method but returns
        structures instead of `MetricsData` objects.
        """
        return [
            data.typed_structure
            for data in self.get_computed_subset(stable, unique, novel, unmatchable)
        ]

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[GenMatStructure]] = OrderedDict(
            [
                ("Stable", self.get_structs_subset(stable=True)),
                ("Unique", self.get_structs_subset(unique=True)),
                ("Novel",  self.get_structs_subset(novel=True)),
                ("Stable Unique", self.get_structs_subset(stable=True, unique=True)),
                ("Stable Novel", self.get_structs_subset(stable=True, novel=True)),
                ("Unique Novel", self.get_structs_subset(unique=True, novel=True)),
                ("S.U.N.", self.get_structs_subset(stable=True, unique=True, novel=True)),
                ("Unmatchable", self.get_structs_subset(unmatchable=True))
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
