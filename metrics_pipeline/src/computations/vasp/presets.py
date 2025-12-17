"""Classes to contain and manage Pymatgen and local VASP configuration presets."""

import os
import typing as tp
from enum import Enum
from dataclasses import dataclass

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import (
    VaspInputSet,
    MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPMetalRelaxSet, MVLRelax52Set, MVLScanRelaxSet,
    MPStaticSet, MPSOCSet, MatPESStaticSet, MPScanStaticSet
)

from metrics_pipeline.src.io import load_yaml_as_dict

# TODO: temporary importation fix, to modify once new packages are built
from metrics_pipeline.src.utils.periodic_table import get_all_valence_electrons


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


@dataclass
class DSolStaticSet(MPRelaxSet):
    """
    Initialize VASP input files for Δ-Sol method computations using
    PBE_54_W_HASH pymatgen set of POTCAR files. Parameters are as 
    described in Δ-Sol method original work by Chan et al. in 2010.
    DFT+U corrections are used as proposed by Jain et al. in 2011.

    References:
        - M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).

        - A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller, 
        K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).

    Args:
        structure (Structure):  The Structure to create inputs for. If None, the input
                                set is initialized without a Structure but one must be
                                set separately before the inputs are generated.

        incar_nelect (float):   The number of electrons to put in the NELECT INCAR tag.
                                In Δ-Sol, several computations with distinct number of 
                                electrons are done, this is a convenient arg to set that.
                                If not given, infers the Δ-Sol N0 electrons calculation 
                                from the given structure.

        **kwargs:               kwargs supported by DictSet.
    
    Raises: ValueError if neither structure nor nelect are given at instanciation time.
    """
    CONFIG = load_yaml_as_dict(
        os.path.join(LOCAL_PRESETS_PATH, "DSolStaticSet.yaml"), on_error="raise"
    )

    def __init__(
            self,
            structure: Structure|None = None,
            incar_nelect: float|None = None,
            **kwargs
        ) -> None:
        """DSolStaticSet init."""
        super().__init__(structure, **kwargs)

        if incar_nelect is None:
            try:
                incar_nelect = get_all_valence_electrons(structure)
            except TypeError as e:
                raise ValueError("Either structure or incar_nelect must be given.") from e

        self.incar_nelect = incar_nelect

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: dict[str, tp.Any] = {"MAGMOM": None, "NELECT": self.incar_nelect}
        return updates
