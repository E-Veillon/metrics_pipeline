#!/usr/bin/python
"""
All specific type aliases used in this library are stored here.
"""

from typing import Union, Sequence, Literal
from pathlib import Path

# PYTHON MATERIALS GENOMICS
from pymatgen.core import Element


PathLike = Union[Path, str]

FormulaLike = Union[str, Sequence[Union[str, int, Element]]]

PMGRelaxSetType = Literal[
    "MITRelaxSet", 
    "MPRelaxSet", 
    "MPScanRelaxSet", 
    "MPHSERelaxSet", 
    "MPMetalRelaxSet", 
    "MVLRelax52Set", 
    "MVLScanRelaxSet"
]

PMGRelaxSet = {
    "MITRelaxSet", 
    "MPRelaxSet", 
    "MPScanRelaxSet", 
    "MPHSERelaxSet", 
    "MPMetalRelaxSet", 
    "MVLRelax52Set", 
    "MVLScanRelaxSet"
}

PMGStaticSetType = Literal[
    "MPStaticSet", 
    "MatPESStaticSet", 
    "MPScanStaticSet"
]

PMGStaticSet = {
    "MPStaticSet", 
    "MatPESStaticSet", 
    "MPScanStaticSet"
}
