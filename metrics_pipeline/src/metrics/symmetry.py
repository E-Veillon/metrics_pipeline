"""Compute Symmetry metric."""

import typing as tp
from collections import OrderedDict

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmetryUndeterminedError

from .metric_base import Metric
from src.utils import raise_or_warn

class Symmetry(Metric):
    """
    Compute Symmetry metric.

    Definition
    ----------
    Structures are considered symmetric if their spacegroup symmetry is not one of the
    triclinic crystal systems, i.e. P1 (no space symmetry element) or P-1 (inversion center only).
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
        Compute Symmetry metric.

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

    def is_symmetric(self, structure: Structure) -> bool:
        """Whether a structure has higher symmetry than triclinic system."""
        try:
            analyzer = SpacegroupAnalyzer(structure, symprec=self.symprec, 
                                         angle_tolerance=self.angleprec)
            is_sym = analyzer.get_crystal_system() != "triclinic"
        except SymmetryUndeterminedError as exc:
            raise_or_warn(self.on_error, type(exc), str(exc))
            return False
        return is_sym

    def _compute(self) -> None:
        if self.workers is not None and self.workers == 0:
            sym_indices = list(map(self.is_symmetric, self.structures))
        else:
            sym_indices = process_map(
                self.is_symmetric,
                self.structures,
                max_workers=self.workers,
                chunksize=min(10, len(self) // 100 + 1),
                desc="Searching structures symmetry"
            )
        self._symmetric_structs = [
            struct for idx, struct in enumerate(self.structures) if sym_indices[idx]
        ]
        self._triclinic_structs = [
            struct for idx, struct in enumerate(self.structures) if not sym_indices[idx]
        ]

    @property
    def symmetric_structs(self) -> list[Structure]:
        """List of structures of higher symmetry than a triclinic system."""
        return self._symmetric_structs

    @property
    def triclinic_structs(self) -> list[Structure]:
        """List of structures of triclinic crystal system."""
        return self._triclinic_structs

    def write_result(self, filename: str, verbose: bool = False) -> None:
        subsets: OrderedDict[str, list[Structure]] = OrderedDict(
            [
                ("Symmetric", self._symmetric_structs),
                ("Triclinic", self._triclinic_structs)
            ]
        )
        self._write_filter_metric_result(filename, subsets, verbose)
