#!/usr/bin/python
"""
Functions that are specific to Δ-Sol method application.

Reference of the method:
- M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
"""

import typing as tp
from enum import Enum
from dataclasses import dataclass

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure

from metrics_pipeline.core.utils.periodic_table import get_all_valence_electrons
from metrics_pipeline.core.utils.visual_iterator import VisualIterator


class DSolCalcType(Enum):
    """
    Enum of Δ-Sol calculation types.
    
    Each run type has an attributed index:
    - Neutral run (nelect) = 0;
    - Main anionic run (nelect + n(best)) = 1;
    - Main cationic run (nelect - n(best)) = 2;
    - Minimally anionic run (nelect + n(min)) = 3;
    - Minimally cationic run (nelect - n(min)) = 4;
    - Maximally anionic run (nelect + n(max)) = 5;
    - Maximally cationic run (nelect - n(max)) = 6;
    """
    NEUTRAL = 0
    BEST_PLUS = 1
    BEST_MINUS = 2
    MIN_PLUS = 3
    MIN_MINUS = 4
    MAX_PLUS = 5
    MAX_MINUS = 6


class DSolFunc(Enum):
    """Enum of Δ-Sol supported DFT functionals."""
    LDA = "LDA"
    PBE = "PBE"
    AM05 = "AM05"


class DSolInput:
    """
    Compute necessary parameters for a given Δ-Sol run on a given structure.

    Reference
    ---------
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
    """
    _EL_PER_XC_VOL: dict[str, dict[str, int]] = {
        "MIN": {
            "LDA_spd": 50, "PBE_spd": 59, "AM05_spd": 60, 
            "LDA_sp": 43, "PBE_sp": 52, "AM05_sp": 52
        }, 
        "BEST": {
            "LDA_spd": 63, "PBE_spd": 72, "AM05_spd": 76, 
            "LDA_sp": 56, "PBE_sp": 68, "AM05_sp": 70
        }, 
        "MAX": {
            "LDA_spd": 80, "PBE_spd": 88, "AM05_spd": 91, 
            "LDA_sp": 78, "PBE_sp": 87, "AM05_sp": 92
        }
    }
    """
    Values of N* (the number of electrons per exchange-correlation volume) 
    used by M.K.Y. Chan and G. Ceder to determine the number "n" of electrons 
    to add to or remove from the simulated structure in the Δ-Sol method.
    See Table I in the original paper for the values.
    """
    def __init__(
        self,
        structure: Structure,
        calc_type: int | DSolCalcType,
        functional: str | DSolFunc
    ) -> None:
        """
        Compute necessary parameters for a given Δ-Sol run on a given structure.

        Parameters
        ----------
        structure: Structure
            Structure to compute band gap of with the Δ-Sol method.

        calc_type: int | DSolCalcType
            Type of run to apply. See `DSolCalcType` class for more details on the
            corresponding indices for each run.

        functional: str | DSolFunc
            One of the DFT functionals supported by the Δ-Sol method. See `DSolFunc` class
            for more details on valid values.
        """
        self.structure = structure
        self.calc_type = DSolCalcType(calc_type)
        self.functional = DSolFunc(functional)

        self.valence_type = self._get_valence_type()
        self.valence_electrons = get_all_valence_electrons(structure)
        self.n_star = self._get_n_star()
        self.n_ratio = self._get_n_ratio()

    def _get_valence_type(self) -> tp.Literal["sp", "spd"]:
        """Get valence orbital types in the structure for Δ-Sol."""
        if self.structure.composition.contains_element_type("f-block"):
            raise NotImplementedError(
                "f-block elements are not supported in Δ-Sol method."
            )
        if self.structure.composition.contains_element_type("d-block"):
            return "spd"

        return "sp"

    def _get_n_star(self) -> int:
        """Get tabulated N* value corresponding to initialized parameters."""
        if self.calc_type == DSolCalcType.NEUTRAL:
            return 0

        value_name = "_".join((self.functional.value, self.valence_type))
        return self._EL_PER_XC_VOL[self.calc_type.name.split(sep="_")[0]][value_name]

    def _get_n_ratio(self) -> float:
        """Get n = N0 / N* corresponding to initialized parameters."""
        if self.calc_type == DSolCalcType.NEUTRAL:
            return 0.0

        return self.valence_electrons / float(self.n_star)

@dataclass
class DSolStructure:
    """
    Dataclass associating a structure with its energy data from Δ-Sol runs.
    
    Attributes
    ----------
    structure: Structure
        Computed structure.

    name: str
        Name of the structure (typically the name of the directory containing all its Δ-Sol runs).

    dft_functional: str | DSolFunc
        Δ-Sol supported functional that was used to do the runs.
        See `DSolFunc` class for more details.

    neutral_energy: float
        Total energy E(nelect) of the neutral structure.

    best_anionic_energy: float
        Total energy E(nelect + n(best)) of the anionized structure
        for band gap determination.

    best_cationic_energy: float
        Total energy E(nelect - n(best)) of the cationized structure
        for band gap determination.

    min_anionic_energy: float
        Total energy E(nelect + n(min)) of the minimally anionized structure
        for uncertainty determination.

    min_cationic_energy: float
        Total energy E(nelect - n(min)) of the minimally cationized structure
        for uncertainty determination.

    max_anionic_energy: float
        Total energy E(nelect + n(max)) of the maximally anionized structure
        for uncertainty determination.

    max_cationic_energy: float
        Total energy E(nelect - n(max)) of the maximally cationized structure
        for uncertainty determination.
    """
    structure: Structure
    name: str
    functional: str
    neutral_energy: float
    best_anionic_energy: float
    best_cationic_energy: float
    min_anionic_energy: float | None = None
    min_cationic_energy: float | None = None
    max_anionic_energy: float | None = None
    max_cationic_energy: float | None = None

    @property
    def band_gap_energies(self) -> dict[str, float]:
        """Dict of energy values used specifically for band gap value computation."""
        return {
            "best_anionic_energy": self.best_anionic_energy,
            "best_cationic_energy": self.best_cationic_energy
        }

    @property
    def uncertainty_energies(self) -> dict[str, float | None]:
        """Dict of energy values used specifically for band gap uncertainty computations."""
        return {
            "min_anionic_energy": self.min_anionic_energy,
            "min_cationic_energy": self.min_cationic_energy,
            "max_anionic_energy": self.max_anionic_energy,
            "max_cationic_energy": self.max_cationic_energy
        }


def get_dsol_band_gap(dsol_struct: DSolStructure) -> tuple[str, dict[str, float]]:
    """
    Calculate Δ-Sol band gap of one structure from computed Δ-Sol energies.

    Parameters
    ----------
    dsol_struct: DSolStructure
        Structure with its extracted energies.

    Returns
    -------
    (str, dict[str, float])
        Name of the structure and dict of its computed band gap values.
    """
    n_ratio_best = DSolInput(
        dsol_struct.structure, DSolCalcType.BEST_PLUS, dsol_struct.functional
    ).n_ratio

    # E_FG = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n -> Δ-Sol band gap
    e_diff_best = (
        dsol_struct.best_anionic_energy
        + dsol_struct.best_cationic_energy
        - 2 * dsol_struct.neutral_energy
    )
    e_bg_best = e_diff_best / n_ratio_best

    if any(energy is None for energy in dsol_struct.uncertainty_energies.values()):
        return dsol_struct.name, {"E_band_gap": e_bg_best}

    # Compute uncertainties
    n_ratio_min = DSolInput(
        dsol_struct.structure, DSolCalcType.MIN_PLUS, dsol_struct.functional
    ).n_ratio
    n_ratio_max = DSolInput(
        dsol_struct.structure, DSolCalcType.MAX_PLUS, dsol_struct.functional
    ).n_ratio
    e_diff_min = (
        tp.cast(float, dsol_struct.min_anionic_energy)
        + tp.cast(float, dsol_struct.min_cationic_energy)
        - 2 * dsol_struct.neutral_energy
    )
    e_bg_min = e_diff_min / n_ratio_min

    e_diff_max = (
        tp.cast(float, dsol_struct.max_anionic_energy)
        + tp.cast(float, dsol_struct.max_cationic_energy)
        - 2 * dsol_struct.neutral_energy
    )
    e_bg_max = e_diff_max / n_ratio_max

    e_band_gap_dict = {
        "E_band_gap": e_bg_best,
        "E_band_gap_min": e_bg_min,
        "E_band_gap_max": e_bg_max
    }

    return dsol_struct.name, e_band_gap_dict


def batch_get_dsol_band_gaps(
        dsol_structs: list[DSolStructure],
        workers: int | None = None
    ) -> dict[str, dict[str, float]]:
    """
    Calculate Δ-Sol band gap of structures from computed Δ-Sol energies.

    Parameters
    ----------
    dsol_structs: list[DSolStructure]
        Structures associated with their extracted energies from Δ-Sol runs.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    Returns
    -------
    dict[str, dict[str, float]]
        Nested dict of {'struct_name': {band gap values}}.
    """
    description = "Computing band gap values"
    if workers == 0:
        e_band_gaps = [
            get_dsol_band_gap(dsol_struct) for dsol_struct in VisualIterator(
                dsol_structs, desc=description, unit="band gaps computed", percent=True
            )
        ]
    else:
        e_band_gaps = process_map(
            get_dsol_band_gap,
            dsol_structs,
            max_workers=workers,
            chunksize=min(10, len(dsol_structs) // 100 + 1),
            desc=description
        )

    return dict(e_band_gaps)
