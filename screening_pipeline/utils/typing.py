'''
All specific type aliases used in this library are stored here.
'''


from typing import Union, Iterable, Literal
from pathlib import Path
from pymatgen.core.structure import Element


PathLike    = Union[Path, str]

FormulaLike = Union[str, Iterable[Union[str, int, Element]]]

PMGRelaxSet = Literal[
    'MITRelaxSet', 
    'MPRelaxSet', 
    'MPScanRelaxSet', 
    'MPHSERelaxSet', 
    'MPMetalRelaxSet', 
    'MVLRelax52Set', 
    'MVLScanRelaxSet'
]

PMGStaticSet = Literal[
    'MPStaticSet', 
    'MatPESStaticSet', 
    'MPScanStaticSet'
]
