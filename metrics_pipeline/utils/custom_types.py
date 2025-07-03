#!/usr/bin/python
"""
All specific type aliases used in this library are stored here.
"""


import typing as tp
from pathlib import Path
from enum import Enum

# PYTHON MATERIALS GENOMICS
from pymatgen.core import Element


PathLike = Path | str

FormulaLike = str | tp.Sequence[str | int | Element]

class PMGRelaxSet(Enum):
    """Enum class of known to date VASP relaxation presets implemented in pymatgen."""
    MITRELAXSET = "MITRelaxSet"
    MPRELAXSET = "MPRelaxSet"
    MPSCANRELAXSET = "MPScanRelaxSet"
    MPMETALRELAXSET = "MPMetalRelaxSet"
    MVLRELAX52SET = "MVLRelax52Set"
    MVLSCANRELAXSET = "MVLScanRelaxSet"

    @property
    def names(self) -> list[str]:
        """Get the list of all Enum members names."""
        return list(member.name for member in self) # type: ignore

    @property
    def values(self) -> list[str]:
        """Get the list of all Enum members values."""
        return list(member.value for member in self) # type: ignore

class PMGStaticSet(Enum):
    """Enum class of known to date VASP static presets implemented in pymatgen."""
    MPSTATICSET = "MPStaticSet"
    MATPESSTATICSET = "MatPESStaticSet"
    MPSCANSTATICSET = "MPScanStaticSet"
    MPSOCSET = "MPSOCSet"

    @property
    def names(self) -> list[str]:
        """Get the list of all Enum members names."""
        return list(member.name for member in self) # type: ignore

    @property
    def values(self) -> list[str]:
        """Get the list of all Enum members values."""
        return list(member.value for member in self) # type: ignore
