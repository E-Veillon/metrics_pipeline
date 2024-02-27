'''
All specific type aliases used in this library are stored here.
'''


from typing import Union, Sequence, Literal
from pathlib import Path
from pymatgen.core.structure import Element
from pymatgen.io.vasp.sets import (
    MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPHSERelaxSet, MPMetalRelaxSet, MVLScanRelaxSet, MVLRelax52Set,
    MPStaticSet, MPScanStaticSet, MatPESStaticSet
)

PathLike = Union[Path, str]

FormulaLike = Union[str, Sequence[Union[str, int, Element]]]

PMGRelaxSetType = Literal[
    'MITRelaxSet', 
    'MPRelaxSet', 
    'MPScanRelaxSet', 
    'MPHSERelaxSet', 
    'MPMetalRelaxSet', 
    'MVLRelax52Set', 
    'MVLScanRelaxSet'
]

PMGRelaxSet = {
    'MITRelaxSet', 
    'MPRelaxSet', 
    'MPScanRelaxSet', 
    'MPHSERelaxSet', 
    'MPMetalRelaxSet', 
    'MVLRelax52Set', 
    'MVLScanRelaxSet'
}

PMGStaticSetType = Literal[
    'MPStaticSet', 
    'MatPESStaticSet', 
    'MPScanStaticSet'
]

PMGStaticSet = {
    'MPStaticSet', 
    'MatPESStaticSet', 
    'MPScanStaticSet'
}
