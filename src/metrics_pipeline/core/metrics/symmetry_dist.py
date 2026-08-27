"""Compute spacegroup symmetry distribution of structures with flexible filtering features."""

import typing as tp
from collections import OrderedDict

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmetryUndeterminedError

from .metric_base import Metric, MetricsData
from core.utils import raise_or_warn, check_type
from core.utils.genmat_data import GenMatStructure
from core.utils.spg_data import (
    ALL_CRYSTAL_FAMILIES, ALL_CRYSTAL_SYSTEMS, ALL_POINT_GROUPS, ALL_SPACEGROUPS, Spacegroup
)
from core.genmat_io import JsonWriter


class SymmetryClassifier(Metric):
    """
    Compute spacegroup symmetry distribution of structures with flexible filtering features. 

    Definition
    ----------
    Spacegroup symmetry of all structures is computed. Filtering methods return subsets
    of structures according to their symmetry, either at spacegroup, point group, crystal system
    or crystal family level.
    """
    _on_error: tp.Literal["raise", "warn", "ignore"]

    def __init__(
        self,
        structures: list[GenMatStructure],
        symprec: float = 0.01,
        angleprec: float = 5.0,
        workers: int | None = None,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> None:
        """
        Compute spacegroup symmetry distribution of structures with flexible filtering features.

        Parameters
        ----------
        structures: list[GenMatStructure]
            Structures to calculate Symmetry on.

        symprec: float
            Tolerance for symmetry finding. See `pymatgen.symmetry.analyzer.SpacegroupAnalyzer`
            for more details. Defaults to 0.01.

        angleprec: float
            Angle tolerance for symmetry finding. See `pymatgen.symmetry.analyzer.SpacegroupAnalyzer`
            for more details. Defaults to 5.0 degrees.

        workers: int, optional
            Number of parallel processes to spawn for high throughput symmetry finding.
            If not given, will use default of `tqdm.contrib.concurrent.process_map()`.
            pass 0 to disable `process_map()` and do it sequentially.

        on_error: "raise", "warn", "ignore"
            What to do in case the symmetry of a structure cannot be found.
            Defaults to "warn".
        """
        super().__init__(structures, workers)
        if on_error not in {"raise", "warn", "ignore"}:
            raise ValueError(
                f"'on_error' only supports 'raise', 'warn' or 'ignore', got {on_error!r}."
        )
        self._symprec = symprec
        self._angleprec = angleprec
        self._on_error = on_error

        self._compute()

    def compute_symmetry(self, structure: GenMatStructure) -> GenMatStructure:
        """
        Compute symmetry of a crystal using pymatgen's SpacegroupAnalyzer
        with initialized parameters.

        Parameters
        ----------
        structure: GenMatStructure
            Structure to compute symmetry for.

        Returns
        -------
        GenMatStructure
            The same structure with symmetry data stored in the `spacegroup` attribute.
        """
        try:
            analyzer = SpacegroupAnalyzer(
                structure, symprec=self._symprec, angle_tolerance=self._angleprec
            )
        except SymmetryUndeterminedError as exc:
            raise_or_warn(self._on_error, type(exc), str(exc))
            spg = Spacegroup(0, is_uncomputable=True)
        else:
            spg = Spacegroup(analyzer.get_space_group_number())
        structure.spacegroup = spg
        return structure

    def get_symmetry_data(self, structure: GenMatStructure) -> MetricsData:
        """
        Compute symmetry of a crystal using pymatgen's SpacegroupAnalyzer
        with initialized parameters. Get the result as a computed MetricsData
        object.

        Parameters
        ----------
        structure: GenMatStructure
            Structure to compute symmetry for.

        Returns
        -------
        MetricsData
            Computed MetricsData containing the analyzed structure and metric
            results.
        """
        sym_struct = self.compute_symmetry(structure)
        return MetricsData(
            sym_struct,
            is_symmetric=sym_struct.spacegroup.int_number > 2,
            additional_data=self._get_metric_settings()
        )

    def _get_metric_settings(self) -> dict[str, tp.Any]:
        """Get a dict of initialized parameters for this metric."""
        return {
            f"{type(self).__name__}_settings": {
                "symprec": self._symprec,
                "angleprec": self._angleprec
            }
        }

    def _compute(self) -> None:
        desc = "Computing structures symmetry"
        if self.workers == 0:
            self._computed_data = list(
                self._sequential_compute(
                    self.get_symmetry_data, self.structures, desc=desc
                )
            )
        else:
            self._computed_data = self._parallel_compute(
                self.get_symmetry_data, self.structures, desc=desc
            )

    @property
    def computed_data(self) -> list[MetricsData]:
        """List of all computed MetricData objects."""
        return self._computed_data

    @property
    def computable_data(self) -> list[MetricsData]:
        """
        List of computed MetricsData containing structures whose spacegroup
        could be determined.
        """
        return [data for data in self.computed_data if not data.typed_structure.spacegroup.is_uncomputable]

    @property
    def uncomputable_data(self) -> list[MetricsData]:
        """
        List of computed MetricsData containing structures whose spacegroup
        could not be determined.
        """
        return [data for data in self.computed_data if data.typed_structure.spacegroup.is_uncomputable]

    @property
    def computable_structs(self) -> list[GenMatStructure]:
        """List of structures whose spacegroup symmetry could be determined."""
        return [data.typed_structure for data in self.computed_data if not data.typed_structure.spacegroup.is_uncomputable]
    
    @property
    def uncomputable_structs(self) -> list[GenMatStructure]:
        """List of structures whose spacegroup symmetry could not be determined."""
        return [data.typed_structure for data in self.computed_data if data.typed_structure.spacegroup.is_uncomputable]

    @property
    def computable_names(self) -> list[str]:
        """List of names of structures whose spacegroup symmetry could be determined."""
        return [data.typed_structure.name for data in self.computed_data if not data.typed_structure.spacegroup.is_uncomputable]

    @property
    def uncomputable_names(self) -> list[str]:
        """List of names of structures whose spacegroup symmetry could not be determined."""
        return [data.typed_structure.name for data in self.computed_data if data.typed_structure.spacegroup.is_uncomputable]

    @property
    def computable_names_set(self) -> set[str]:
        """Set of names of structures whose spacegroup symmetry could be determined."""
        return {data.typed_structure.name for data in self.computed_data if not data.typed_structure.spacegroup.is_uncomputable}

    @property
    def uncomputable_names_set(self) -> set[str]:
        """Set of names of structures whose spacegroup symmetry could not be determined."""
        return {data.typed_structure.name for data in self.computed_data if data.typed_structure.spacegroup.is_uncomputable}

    def get_symmetry_subset(
        self, symmetry: str | int, family: bool = False
    ) -> list[GenMatStructure]:
        """
        Get the subset of structures of a particular symmetry class.

        Parameters
        ----------
        symmetry: str | int
            Name or symbol of the crystal family, systeem, point group or space group,
            or integer index of the space group to filter by (e.g. 'cubic', 'm-3m',
            'Fm-3m', or 225 are all valid inputs).

        family: bool
            Whether passed crystal class name refers to the crystal system or crystal family.
            Only used if 'hexagonal' is passed to avoid ambiguity between the hexagonal system
            and the hexagonal family, which includes both hexagonal and trigonal systems.
            Defaults to False.

        Raises
        ------
        `ValueError`: If `symmetry` is not a valid symmetry class.

        Returns
        -------
        list[GenMatStructure]
            List of structures of given symmetry.
        """
        check_type(symmetry, "symmetry", (str, int))

        if isinstance(symmetry, int) and 1 <= symmetry <= 230:
            return [
                struct for struct in self.computable_structs
                if struct.spacegroup.int_number == symmetry
            ]
        if symmetry in ALL_SPACEGROUPS:
            return [
                struct for struct in self.computable_structs
                if struct.spacegroup.symbol == symmetry
            ]
        if symmetry in ALL_POINT_GROUPS:
            return [
                struct for struct in self.computable_structs
                if struct.spacegroup.point_group == symmetry
            ]
        if not family and symmetry in ALL_CRYSTAL_SYSTEMS:
            return [
                struct for struct in self.computable_structs
                if struct.spacegroup.crystal_system == symmetry
            ]
        if family and symmetry in ALL_CRYSTAL_FAMILIES:
            return [
                struct for struct in self.computable_structs
                if struct.spacegroup.crystal_family == symmetry
            ]
        if family and symmetry == "trigonal":
            raise ValueError(
                f"The {symmetry!r} system was passed but 'family' was set to True. "
                f"Set 'family' to False to get the {symmetry} system subset."
            )
        raise ValueError(
            f"{symmetry!r} is not a valid symmetry class name or symbol or space group index."
        )

    def write_result(self, filename: str, verbose: bool = False) -> None:
        system_subsets = OrderedDict(
            [
                ("Uncomputable", self.uncomputable_structs),
                ("Triclinic", self.get_symmetry_subset("triclinic")),
                ("Monoclinic", self.get_symmetry_subset("monoclinic")),
                ("Orthorhombic", self.get_symmetry_subset("orthorhombic")),
                ("Tetragonal", self.get_symmetry_subset("tetragonal")),
                ("Trigonal", self.get_symmetry_subset("trigonal")),
                ("Hexagonal", self.get_symmetry_subset("hexagonal")),
                ("Cubic", self.get_symmetry_subset("cubic"))
            ]
        )
        pg_subsets = OrderedDict(
            [
                (f"Point Group {pg!r}", self.get_symmetry_subset(pg))
                for pg in ALL_POINT_GROUPS
            ]
        )
        spg_subsets = OrderedDict(
            [
                (f"Spacegroup {spg!r} ({num})", self.get_symmetry_subset(num))
                for num, spg in enumerate(ALL_SPACEGROUPS, start=1)
            ]
        )
        subsets = system_subsets | pg_subsets | spg_subsets
        self._write_filter_metric_result(filename, subsets, verbose)

    def write_json(self, filename: str, verbose: bool = False) -> None:
        """
        Write a compact JSON file containing all results of this metric
        at crystal system, point group and space group levels.

        Parameters
        ----------
        filename: str
            Path to output JSON file to write.

        verbose: bool
            Whether to add the list of headers in each symmetry class.
            Defaults to False.
        """
        results = {
            "uncomputable": {}, "by_system": {}, "by_point_group": {}, "by_space_group": {}
        }
        results["uncomputable"] = {"total": len(self.uncomputable_structs)}
        if verbose:
            results["uncomputable"] |= {
                "headers": [struct.name for struct in self.uncomputable_structs]
            }
        for system in ALL_CRYSTAL_SYSTEMS:
            structures = self.get_symmetry_subset(system)
            results["by_system"][system] = {"total": len(structures)}
            if verbose:
                results["by_system"][system] |= {
                    "headers": [struct.name for struct in structures]
                }
        for pg in ALL_POINT_GROUPS:
            structures = self.get_symmetry_subset(pg)
            results["by_point_group"][system] = {"total": len(structures)}
            if verbose:
                results["by_point_group"][system] |= {
                    "headers": [struct.name for struct in structures]
                }
        for sg in ALL_SPACEGROUPS:
            structures = self.get_symmetry_subset(sg)
            results["by_space_group"][system] = {"total": len(structures)}
            if verbose:
                results["by_space_group"][system] |= {
                    "headers": [struct.name for struct in structures]
                }
        JsonWriter(filename, results, indent=4).write_as_dict()
