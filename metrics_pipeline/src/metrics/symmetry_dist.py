"""Compute spacegroup symmetry distribution of structures with flexible filtering features."""

import typing as tp
from collections import OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmetryUndeterminedError

from .metric_base import Metric
from src.utils import raise_or_warn
from src.utils.spg_data import (
    ALL_CRYSTAL_FAMILIES, ALL_CRYSTAL_SYSTEMS, ALL_POINT_GROUPS, ALL_SPACEGROUPS
)
from src.io import JsonWriter

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
            under the key defined by the `sym_data_key` property. The property contains
            `None` instead of a dict if symmetry could not be determined and `on_error`
            was not set to `raise`.
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
                self.sym_data_sg_idx_key: analyzer.get_space_group_number(),
                self.sym_data_sg_key: analyzer.get_space_group_symbol(),
                self.sym_data_pg_key: analyzer.get_point_group_symbol(),
                self.sym_data_system_key: cs,
                self.sym_data_family_key: "hexagonal" if cs == "trigonal" else cs
            }

        return structure

    def _compute(self) -> None:
        if self.workers == 0:
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
    def all_sym_classes(self) -> dict[str, tuple[str, ...]]:
        """
        All valid symmetry classes in a dict of {'symmetry class': (all symbols in the class)}.
        """
        return {
            "family": ALL_CRYSTAL_FAMILIES,
            "system": ALL_CRYSTAL_SYSTEMS,
            "point_group": ALL_POINT_GROUPS,
            "space_group": ALL_SPACEGROUPS
        }

    @property
    def sym_data_key(self) -> str:
        """Key in `properties` attribute of structures where symmetry data is stored."""
        return "GenMat_symmetry_data"

    @property
    def sym_data_family_key(self) -> str:
        return "crystal_family"
    
    @property
    def sym_data_system_key(self) -> str:
        return "crystal_system"

    @property
    def sym_data_pg_key(self) -> str:
        return "point_group_symbol"

    @property
    def sym_data_sg_key(self) -> str:
        return "space_group_symbol"

    @property
    def sym_data_sg_idx_key(self) -> str:
        return "space_group_number"

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
        structures properties is set to `None`.
        """
        return self._uncomputable_structs

    def get_symmetry_subset(
        self, sym_class: str, symmetry: str | int, strict: bool = True
    ) -> list[Structure]:
        """
        Get structures subset of a particular symmetry.

        Parameters
        ----------
        sym_class: str
            The symmetry class to filter by. Can be either 'family', 'system', 'point_group'
            or 'space_group'. Notably used to avoid ambiguities such as between the hexagonal
            crystal system and the hexagonal family which contains both hexagonal and trigonal
            systems.

        symmetry: str | int
            Name or symbol of the symmetry, or integer index of the space group to filter by
            (e.g. 'cubic', 'm-3m', 'Fm-3m', or 225 are all valid inputs).

        strict: bool
            Whether to raise a `ValueError` (True) or return an empty list (False) if one of
            passed arguments is invalid. Defaults to True.

        Raises
        ------
        `ValueError`: If one of `sym_class` or `symmetry` passed arguments are invalid and `strict`
        is set to `True`.

        Returns
        -------
        list[Structure]
            List of structures of given symmetry.
        """
        # Define internal error checker
        def _raise_or_empty(msg: str) -> list:
            if strict:
                raise ValueError(msg)
            return []

        # Deal with integer case
        if isinstance(symmetry, int):
            if symmetry <= 0 or symmetry > 230:
                return _raise_or_empty(
                    f"Space group indices are between 1 and 230, got {symmetry}."
                )
            return [
                    struct for struct in self._computed_structs
                    if struct.properties[self.sym_data_key][self.sym_data_sg_idx_key] == symmetry
                ]

        # Check validity of symmetry arguments
        if not symmetry in sum((ALL_CRYSTAL_SYSTEMS, ALL_POINT_GROUPS, ALL_SPACEGROUPS), tuple()):
            return _raise_or_empty(
                    f"{symmetry!r} is not a recognized symmetry class name or symbol."
                )
        sym_types = self.all_sym_classes.get(sym_class)
        if sym_types is None:
            sym_classes = ", ".join(list(self.all_sym_classes.keys()))
            return _raise_or_empty(
                f"{sym_class!r} si not a valid symmetry class, must be either of {sym_classes}."
            )
        if not symmetry in sym_types:
            class_name = (
                f"crystal {sym_class}" if sym_class in {"family, system"}
                else f"{sym_class.replace('_', '')} symbol"
            )
            return _raise_or_empty(f"{symmetry!r} is not a valid {class_name}.")

        # Filter structures with defined symmetry
        sym_class_to_key = {
            "family": self.sym_data_family_key,
            "system": self.sym_data_system_key,
            "point_group": self.sym_data_pg_key,
            "space_group": self.sym_data_sg_key
        }
        return [
            struct for struct in self._computed_structs
            if struct.properties[self.sym_data_key][sym_class_to_key[sym_class]] == symmetry
        ]

    def write_result(self, filename: str, verbose: bool = False) -> None:
        system_subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Uncomputable", self.uncomputable_structs),
                ("Triclinic", self.get_symmetry_subset("system", "triclinic")),
                ("Monoclinic", self.get_symmetry_subset("system", "monoclinic")),
                ("Orthorhombic", self.get_symmetry_subset("system", "orthorhombic")),
                ("Tetragonal", self.get_symmetry_subset("system", "tetragonal")),
                ("Trigonal", self.get_symmetry_subset("system", "trigonal")),
                ("Hexagonal", self.get_symmetry_subset("system", "hexagonal")),
                ("Cubic", self.get_symmetry_subset("system", "cubic"))
            ]
        )
        pg_subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                (f"Point Group {pg!r}", self.get_symmetry_subset("point_group", pg))
                for pg in ALL_POINT_GROUPS
            ]
        )
        spg_subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                (f"Spacegroup {spg!r} ({num})", self.get_symmetry_subset("space_group", num))
                for num, spg in enumerate(ALL_SPACEGROUPS, start=1)
            ]
        )
        subsets: OrderedDict[str, list[Structure]] = system_subsets | pg_subsets | spg_subsets
        self._write_filter_metric_result(filename, subsets, verbose)

    def write_json(self, filename: str, verbose: bool = False) -> None:
        """
        Write compact JSON file containing all results of this metric
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

        results["uncomputable"] = {
                "total": len(self.uncomputable_structs),
                "headers": [struct.properties["header"] for struct in self.uncomputable_structs]
            } if verbose else {"total": len(self.uncomputable_structs)}

        for system in ALL_CRYSTAL_SYSTEMS:
            structures = self.get_symmetry_subset("system", system)
            results["by_system"][system] = {
                "total": len(structures),
                "headers": [struct.properties["header"] for struct in structures]
            } if verbose else {"total": len(structures)}

        for pg in ALL_POINT_GROUPS:
            structures = self.get_symmetry_subset("point_group", pg)
            results["by_point_group"][pg] = {
                "total": len(structures),
                "headers": [struct.properties["header"] for struct in structures]
            } if verbose else {"total": len(structures)}

        for sg in ALL_SPACEGROUPS:
            structures = self.get_symmetry_subset("space_group", sg)
            results["by_space_group"][sg] = {
                "total": len(structures),
                "headers": [struct.properties["header"] for struct in structures]
            } if verbose else {"total": len(structures)}

        JsonWriter(filename, results, indent=4).write_as_dict()
