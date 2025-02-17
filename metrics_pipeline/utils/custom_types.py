#!/usr/bin/python
"""
All specific type aliases used in this library are stored here.
"""

from typing import Union, Sequence, List
from pathlib import Path
from enum import Enum

# PYTHON MATERIALS GENOMICS
from pymatgen.core import Element


PathLike = Union[Path, str]

FormulaLike = Union[str, Sequence[Union[str, int, Element]]]

class PMGRelaxSet(Enum):
    """Enum class of known to date VASP relaxation presets implemented in pymatgen."""
    MITRELAXSET = "MITRelaxSet"
    MPRELAXSET = "MPRelaxSet"
    MPSCANRELAXSET = "MPScanRelaxSet"
    MPMETALRELAXSET = "MPMetalRelaxSet"
    MVLRELAX52SET = "MVLRelax52Set"
    MVLSCANRELAXSET = "MVLScanRelaxSet"

    @property
    def names(self) -> List[str]:
        """Get the list of all Enum members names."""
        return list(member.name for member in self)

    @property
    def values(self) -> List[str]:
        """Get the list of all Enum members values."""
        return list(member.value for member in self)

class PMGStaticSet(Enum):
    """Enum class of known to date VASP static presets implemented in pymatgen."""
    MPSTATICSET = "MPStaticSet"
    MATPESSTATICSET = "MatPESStaticSet"
    MPSCANSTATICSET = "MPScanStaticSet"
    MPSOCSET = "MPSOCSet"

    @property
    def names(self) -> List[str]:
        """Get the list of all Enum members names."""
        return list(member.name for member in self)

    @property
    def values(self) -> List[str]:
        """Get the list of all Enum members values."""
        return list(member.value for member in self)
