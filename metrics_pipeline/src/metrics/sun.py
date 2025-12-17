"""Compute the Stability, Unicity, Novelty compound metric (S.U.N.)."""

from collections import OrderedDict

from pymatgen.core import Structure

from .metric_base import Metric
from .stability import Stability
from .unicity import Unicity
from .novelty import Novelty


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
    def __init__(
        self, structures: list[Structure], ref_structs: list[Structure], **kwargs
    ) -> None:
        """
        Compute the Stability, Unicity, Novelty compound metric (S.U.N.).
        
        Parameters
        ----------
        structures: list[Structure]
            Structures to calculate S.U.N. metrics on.

        ref_structs: list[Structure]
            Structures of known formation energy to use as references.

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
        """
        super().__init__(structures)
        self.ref_structs = ref_structs
        self.stable_tol = kwargs.pop("stable_tol", 0.1)
        self.ltol = kwargs.pop("ltol", 0.2)
        self.stol = kwargs.pop("stol", 0.3)
        self.angle_tol = kwargs.pop("angle_tol", 5.0)

        self._compute()

    def _compute(self) -> None:
        stability = Stability(self.structures, self.ref_structs, self.stable_tol)
        unicity = Unicity(self.structures, self.ltol, self.stol, self.angle_tol)
        novelty = Novelty(self.structures, self.ref_structs, self.ltol, self.stol, self.angle_tol)
        for struct in self.structures:
            struct.properties[self.stable_key] = struct in stability.stable_structs
            struct.properties[self.unique_key] = (
                struct in unicity.unique_structs + unicity.unmatchable_structs
            )
            struct.properties[self.novel_key] = (
                struct in novelty.novel_structs + novelty.unmatchable_structs
            )
            struct.properties[self.unmatch_key] = (
                struct in unicity.unmatchable_structs + novelty.unmatchable_structs
            )
        self._computed_structs = self.structures

    @property
    def stable_key(self) -> str:
        """
        Key in the computed structures properties where their Stability status is stored.
        """
        return "SUN_is_stable"

    @property
    def unique_key(self) -> str:
        """
        Key in the computed structures properties where their Unicity status is stored.
        """
        return "SUN_is_unique"

    @property
    def novel_key(self) -> str:
        """
        Key in the computed structures properties where their Novelty status is stored.
        """
        return "SUN_is_novel"

    @property
    def unmatch_key(self) -> str:
        """
        Key in the computed structures properties where their unmatchability status is stored.
        """
        return "SUN_is_unmatchable"

    @property
    def computed_structs(self) -> list[Structure]:
        """
        Full list of computed structures with their properties populated with metrics results.
        """
        return self._computed_structs

    def get_computed_subset(
        self,
        stable: bool | None = True,
        unique: bool | None = True,
        novel: bool | None = True,
        unmatchable: bool | None = False
    ) -> list[Structure]:
        """
        Get a subset of computed structures according to given metrics flags.

        Parameters
        ----------
        stable: bool or None
            - If set to True, only return structures that passed the Stability metric.
            - If set to False, only return structures that failed the Stability metric.
            - If set to None, return structures indifferently for this metric.
            Defaults to True.

        unique: bool or None
            - If set to True, only return structures that passed the Unicity metric.
            - If set to False, only return structures that failed the Unicity metric.
            - If set to None, return structures indifferently for this metric.
            Defaults to True.

        novel: bool or None
            - If set to True, only return structures that passed the Novelty metric.
            - If set to False, only return structures that failed the Novelty metric.
            - If set to None, return structures indifferently for this metric.
            Defaults to True.

        unmatchable: bool or None
            - If set to True, only return structures that could not be matched.
            - If set to False, only return structures that could be matched.
            - If set to None, return structures indifferently for this metric.
            Defaults to False.

        Returns
        -------
        list[Structure]
            The list of computed structures with the queried metrics properties.

        Notes
        -----
        - The default behavior returns matched S.U.N. structures only.
        - Passing `None` to all metrics is equivalent to the `computed_structs`
        property return.
        """
        if stable is None and unique is None and novel is None and unmatchable is None:
            return self.computed_structs

        return list(
            filter(
                lambda s: (
                    (True if stable is None else s.properties[self.stable_key] is stable) and
                    (True if unique is None else s.properties[self.unique_key] is unique) and
                    (True if novel is None else s.properties[self.novel_key] is novel) and
                    (True if unmatchable is None else s.properties[self.unmatch_key] is unmatchable)
                ), self.computed_structs
            )
        )

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Stable", self.get_computed_subset(unique=None, novel=None)),
                ("Unique", self.get_computed_subset(stable=None, novel=None)),
                ("Novel",  self.get_computed_subset(stable=None, unique=None)),
                ("Stable Unique", self.get_computed_subset(novel=None)),
                ("Stable Novel", self.get_computed_subset(unique=None)),
                ("Unique Novel", self.get_computed_subset(stable=None)),
                ("S.U.N.", self.get_computed_subset()),
                ("Unmatchable", self.get_computed_subset(
                    stable=None, unique=None, novel=None, unmatchable=True))
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
