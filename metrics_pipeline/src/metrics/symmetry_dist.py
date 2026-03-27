"""Compute spacegroup symmetry distribution of structures with flexible filtering features."""

import typing as tp
from collections import OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmetryUndeterminedError

from .metric_base import Metric
from src.utils import raise_or_warn, ALL_POINT_GROUPS, ALL_SPACEGROUPS

class SymmetryClassifier(Metric):
    """
    Compute spacegroup symmetry distribution of structures with flexible filtering features. 

    Definition
    ----------
    Spacegroup symmetry of all structures is computed. Filtering methods return subsets
    of structures according to their symmetry, either at spacegroup, point group, crystal system
    or crystal family level.
    """
    on_error: tp.Literal["raise", "warn", "ignore"]

    def __init__(
        self,
        structures: list[Structure],
        symprec: float = 0.01,
        angleprec: float = 5.0,
        workers: int | None = None,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> None:
        """
        Compute spacegroup symmetry distribution of structures with flexible filtering features.

        Parameters
        ----------
        structures: list[Structure]
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
        super().__init__(structures)
        if workers is not None:
            if not isinstance(workers, int):
                raise TypeError(
                    f"'workers' expected a type 'int', got  {type(workers).__name__}."
            )
            if workers < 0:
                raise ValueError(
                    f"'workers' must be positive or zero, got {workers}."
            )
        if on_error.lower() not in {"raise", "warn", "ignore"}:
            raise ValueError(
                f"'on_error' only supports 'raise', 'warn' or 'ignore', got {on_error}."
        )
        self.symprec = symprec
        self.angleprec = angleprec
        self.workers = workers
        self.on_error = on_error

        self._compute()

    def get_symmetry(self, structure: Structure) -> Structure:
        """
        Compute symmetry of a crystal using pymatgen's SpacegroupAnalyzer
        with initialized parameters.

        Parameters
        ----------
        structure: Structure
            Structure to compute symmetry for.

        Returns
        -------
        Structure
            The same structure with symmetry data stored in `properties` attribute
            under the key defined by the `sym_data_key` class attribute.
        """
        try:
            analyzer = SpacegroupAnalyzer(structure, symprec=self.symprec, 
                                         angle_tolerance=self.angleprec)
        except SymmetryUndeterminedError as exc:
            raise_or_warn(self.on_error, type(exc), str(exc))
            structure.properties[self.sym_data_key] = None
        else:
            cs = analyzer.get_crystal_system()
            structure.properties[self.sym_data_key] = {
                "space_group_number": analyzer.get_space_group_number(),
                "space_group_symbol": analyzer.get_space_group_symbol(),
                "point_group_symbol": analyzer.get_point_group_symbol(),
                "crystal_system": cs,
                "crystal_family": "hexagonal" if cs == "trigonal" else cs
            }

        return structure

    def _compute(self) -> None:
        if self.workers is not None and self.workers == 0:
            sym_structs = list(map(self.get_symmetry, self.structures))
        else:
            sym_structs = process_map(
                self.get_symmetry,
                self.structures,
                max_workers=self.workers,
                chunksize=min(10, len(self) // 100 + 1),
                desc="Searching structures symmetry"
            )
        self._computed_structs = [
            struct for struct in sym_structs if struct.properties[self.sym_data_key] is not None
        ]
        self._uncomputable_structs = [
            struct for struct in sym_structs if struct.properties[self.sym_data_key] is None
        ]

    @property
    def sym_data_key(self) -> str:
        """Key in `properties` attribute of structures where symmetry data is stored."""
        return "GenMat_symmetry_data"

    @property
    def computed_structs(self) -> list[Structure]:
        """
        List of successfully processed structures.
        Symmetry data is stored in `properties` attribute of each structure,
        under the key defined by the `sym_data_key` class attribute.
        """
        return self._computed_structs

    @property
    def uncomputable_structs(self) -> list[Structure]:
        """
        List of structures for which symmetry could not be computed.
        The key defined by the `sym_data_key` class attribute in these
        structures properties is set to None.
        """
        return self._uncomputable_structs

    def get_crystal_family(self, family: str) -> list[Structure]:
        """
        Filter structures by their crystal family.

        Parameters
        ----------
        family: str
            The crystal family to filter by.

        Returns
        -------
        list[Structure]
            List of structures that satisfy the filtering criteria defined by `filter_func`.
        """
        return [
            struct for struct in self._computed_structs
            if struct.properties[self.sym_data_key]["crystal_family"] == family
        ]

    def get_crystal_system(self, system: str) -> list[Structure]:
        """
        Filter structures by their crystal system.

        Parameters
        ----------
        system: str
            The crystal system to filter by.

        Returns
        -------
        list[Structure]
            List of structures that satisfy the filtering criteria defined by `filter_func`.
        """
        return [
            struct for struct in self._computed_structs
            if struct.properties[self.sym_data_key]["crystal_system"] == system
        ]

    def get_point_group(self, point_group: str) -> list[Structure]:
        """
        Filter structures by their point group.

        Parameters
        ----------
        point_group: str
            The point group to filter by.

        Returns
        -------
        list[Structure]
            List of structures that satisfy the filtering criteria defined by `filter_func`.
        """
        return [
            struct for struct in self._computed_structs
            if struct.properties[self.sym_data_key]["point_group_symbol"] == point_group
        ]

    def get_space_group(self, space_group: int | str) -> list[Structure]:
        """
        Filter structures by their space group.

        Parameters
        ----------
        space_group: str
            The space group to filter by.

        Returns
        -------
        list[Structure]
            List of structures that satisfy the filtering criteria defined by `filter_func`.
        """
        return [
            struct for struct in self._computed_structs
            if space_group in {struct.properties[self.sym_data_key]["space_group_number"], 
                               struct.properties[self.sym_data_key]["space_group_symbol"]}
        ]

    def write_result(self, filename: str, verbose: bool = False) -> None:
        system_subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Uncomputable", self.uncomputable_structs),
                ("Triclinic", self.get_crystal_system("triclinic")),
                ("Monoclinic", self.get_crystal_system("monoclinic")),
                ("Orthorhombic", self.get_crystal_system("orthorhombic")),
                ("Tetragonal", self.get_crystal_system("tetragonal")),
                ("Trigonal", self.get_crystal_system("trigonal")),
                ("Hexagonal", self.get_crystal_system("hexagonal")),
                ("Cubic", self.get_crystal_system("cubic"))
            ]
        )
        pg_subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                (f"Point Group {pg!r}", self.get_point_group(pg)) for pg in ALL_POINT_GROUPS
            ]
        )
        spg_subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                (f"Spacegroup {spg!r} ({num})", self.get_space_group(num))
                for num, spg in enumerate(ALL_SPACEGROUPS, start=1)
            ]
        )
        subsets: OrderedDict[str, list[Structure]] = system_subsets | pg_subsets | spg_subsets
        self._write_filter_metric_result(filename, subsets, verbose)
