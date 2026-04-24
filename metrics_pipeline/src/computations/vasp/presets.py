"""Classes to contain and manage Pymatgen and local VASP configuration presets."""

import os
import typing as tp
from enum import Enum
from dataclasses import dataclass, field

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import (
    VaspInputSet,
    MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPMetalRelaxSet, MVLRelax52Set, MVLScanRelaxSet,
    MPStaticSet, MPSOCSet, MatPESStaticSet, MPScanStaticSet
)

from src.io import load_yaml_as_dict

# TODO: temporary importation fix, to modify once new packages are built
from src.utils.periodic_table import get_all_valence_electrons


LOCAL_PRESETS_PATH = os.path.join(os.path.dirname(__file__), "presets_data")


class PresetEnum(Enum):
    """Base Enum with convenient methods to navigate between stored presets."""
    @classmethod
    def names(cls) -> list[str]:
        """List of all presets names."""
        return list(member.name for member in cls)

    @classmethod
    def values(cls) -> list[type[VaspInputSet]]:
        """List of all preset classes."""
        return list(member.value for member in cls)

    @classmethod
    def items(cls) -> list[tuple[str, type[VaspInputSet]]]:
        """List of name / preset pairs as tuples, similar to the dict.items() method."""
        return [(member.name, member.value) for member in cls]

    @classmethod
    def is_preset(cls, name: str) -> bool:
        """Whether given name corresponds to an existing preset (case insensitive)."""
        return name.upper() in [preset_name.upper() for preset_name in cls.names()]

    @classmethod
    def get_preset(cls, name: str) -> type[VaspInputSet]:
        """Get one of the stored callable preset classes by its name (case insensitive)."""
        if not cls.is_preset(name):
            raise ValueError(
                f"{name!r} is not a valid preset of {cls.__name__}. "
                f"Available presets in this class: {', '.join(cls.names())}."
            )
        for preset_name, preset in cls.items():
            if name.upper() == preset_name.upper():
                return preset

        raise RuntimeError(
            "Preset name was recognized by 'is_preset' but not found when looping. "
            "This error should never occur and is likely a bug. Please open an issue "
            "on the development repo of the library about this."
        )


class PMGRelaxSet(PresetEnum):
    """Enum class of known to date VASP relaxation presets implemented in pymatgen."""
    MITRELAXSET = MITRelaxSet
    MPRELAXSET = MPRelaxSet
    MPSCANRELAXSET = MPScanRelaxSet
    MPMETALRELAXSET = MPMetalRelaxSet
    MVLRELAX52SET = MVLRelax52Set
    MVLSCANRELAXSET = MVLScanRelaxSet


class PMGStaticSet(PresetEnum):
    """Enum class of known to date VASP static presets implemented in pymatgen."""
    MPSTATICSET = MPStaticSet
    MATPESSTATICSET = MatPESStaticSet
    MPSCANSTATICSET = MPScanStaticSet
    MPSOCSET = MPSOCSet


def _get_genmat_settings() -> dict[str, tp.Any]:
    """Load GenMat settings to override pymatgen defaults."""
    return load_yaml_as_dict(
        os.path.join(LOCAL_PRESETS_PATH, "default_settings.yaml"), on_error="raise"
    )


@dataclass
class GenMatRelaxSet(MPRelaxSet):
    """
    VASP input set derived from pymatgen MPRelaxSet for GenMat relaxation step.
    """
    def __post_init__(self) -> None:
        self.CONFIG: dict[str, tp.Any] = super().CONFIG
        self._GENMAT_SETTINGS: dict[str, tp.Any] = _get_genmat_settings()

        for key in self.CONFIG:
            if isinstance(self.CONFIG[key], dict):
                self.CONFIG[key].update(self._GENMAT_SETTINGS.get(key, {}))
            elif key in self._GENMAT_SETTINGS:
                self.CONFIG[key] = self._GENMAT_SETTINGS[key]

        super().__post_init__()

    @property
    def incar_updates(self) -> dict[str, tp.Any]:
        """Get updates to the INCAR config for this calculation type."""
        updates = super().incar_updates
        local_incar: dict[str, tp.Any] = self._GENMAT_SETTINGS.get("INCAR", {})
        for key in {"MAGMOM", "LDAUU", "LDAUL", "LDAUJ"}:
            local_incar.pop(key, None)
        updates.update(local_incar)
        return updates


class GenMatStaticSet(MPStaticSet):
    """
    VASP input set derived from pymatgen MPStaticSet for GenMat static step.
    """
    def __post_init__(self) -> None:
        self.CONFIG: dict[str, tp.Any] = super().CONFIG
        self._GENMAT_SETTINGS: dict[str, tp.Any] = _get_genmat_settings()

        for key in self.CONFIG:
            if isinstance(self.CONFIG[key], dict):
                self.CONFIG[key].update(self._GENMAT_SETTINGS.get(key, {}))
            elif key in self._GENMAT_SETTINGS:
                self.CONFIG[key] = self._GENMAT_SETTINGS[key]

        super().__post_init__()

    @property
    def incar_updates(self) -> dict[str, tp.Any]:
        """Get updates to the INCAR config for this calculation type."""
        updates = super().incar_updates
        local_incar: dict[str, tp.Any] = self._GENMAT_SETTINGS.get("INCAR", {})
        for key in {"MAGMOM", "LDAUU", "LDAUL", "LDAUJ"}:
            local_incar.pop(key, None)
        updates.update(local_incar)
        return updates


@dataclass
class DSolStaticSet(MPStaticSet):
    """
    Initialize VASP input files for Δ-Sol method computations using
    PBE_54_W_HASH pymatgen set of POTCAR files. Parameters are as 
    described in Δ-Sol method original work by Chan et al. in 2010.
    DFT+U corrections are used as proposed by Jain et al. in 2011.

    References
    ----------
    Δ-Sol method:
    - M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).

    DFT + U corrections:
    - A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller, 
    K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).

    Parameters
    ----------
    structure: Structure, optional
        The Structure to create inputs for. If `None`, the input set is initialized without
        a Structure but one must be set separately before the inputs are generated.

    incar_nelect: float, optional
        The number of electrons to put in the `NELECT` INCAR tag. In Δ-Sol, several computations
        with distinct number of electrons are done, this is a convenient arg to set that.
        If not given, infers the Δ-Sol N0 electrons calculation from the given structure.

    kwargs: Any
        keyword arguments supported by `VaspInputSet`.
    
    Raises
    ------
    `ValueError` if neither structure nor nelect are given at initialization time.
    """
    incar_nelect: float | None = None
    CONFIG = load_yaml_as_dict(
        os.path.join(LOCAL_PRESETS_PATH, "DSolStaticSet.yaml"), on_error="raise"
    )

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: dict[str, tp.Any] = super().incar_updates
        if self.incar_nelect is None:
            try:
                self.incar_nelect = self.nelect
            except RuntimeError as e:
                raise ValueError(
                    "Either 'structure' or 'incar_nelect' must be given at initialization."
                ) from e

        updates.update({"MAGMOM": None, "NELECT": self.incar_nelect})
        return updates


class GenMatSet(PresetEnum):
    """Enum class of VASP presets defined for the GenMat-metrics project."""
    GENMATRELAXSET = GenMatRelaxSet
    GENMATSTATICSET = GenMatStaticSet
    DSOLSTATICSET = DSolStaticSet

# Shortcut constants for easy checks and messages
ALL_PRESETS_NAMES = PMGStaticSet.names() + PMGRelaxSet.names() + GenMatSet.names()
ALL_RELAX_PRESETS_NAMES = PMGRelaxSet.names() + [GenMatSet.GENMATRELAXSET.name]
ALL_STATIC_PRESETS_NAMES = PMGStaticSet.names() + [
    GenMatSet.GENMATSTATICSET.name, GenMatSet.DSOLSTATICSET.name
]
STATIC_PRESETS_NAMES_NO_DSOL = PMGStaticSet.names() + [GenMatSet.GENMATSTATICSET.name]

ALL_PRESETS_NAMES_LOWER = [name.lower() for name in ALL_PRESETS_NAMES]
ALL_RELAX_PRESETS_NAMES_LOWER = [name.lower() for name in ALL_RELAX_PRESETS_NAMES]
ALL_STATIC_PRESETS_NAMES_LOWER = [name.lower() for name in ALL_STATIC_PRESETS_NAMES]
STATIC_PRESETS_NAMES_NO_DSOL_LOWER = [name.lower() for name in STATIC_PRESETS_NAMES_NO_DSOL]


if __name__ == "__main__":
    # Test for the DSolStaticSet class. You may have to change the given path
    # to one pointing at a valid CIF structure file for it to work properly.
    # You'll also need to set PMG_VASP_PSP_DIR for POTCAR files in .pmgrc.yaml.
    file_path = ""
    PATHTEST = os.path.join(os.path.expanduser("~"), file_path)
    with open(PATHTEST, "rt", encoding="utf-8") as test_file:
        struct = CifParser(test_file).parse_structures()[0] # type: ignore
    dset = DSolStaticSet(struct).get_input_set()
    with open(
        os.path.join(os.path.dirname(PATHTEST), "DeltaVaspInput.txt"),
        mode="wt", encoding="utf-8"
    ) as out:
        out.write(str(dset))
    print(dset.incar_nelect) # type: ignore
    dset.incar_nelect = 58 # type: ignore
    print(dset.incar_nelect) # type: ignore
