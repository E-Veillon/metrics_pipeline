#!/usr/bin/python
"""
All specific type aliases used in this library are stored here.
"""

import typing as tp
import typing_extensions as tpe
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

class VisualIterator:
    """
    A minimalist way to visualize progression while iterating through long sequences.
    Possibility to custom printed messages.
    """
    def __init__(
        self: tpe.Self,
        iterable: tp.Sequence,
        desc: str | None = None,
        end_desc: str | None = None
    ) -> None:
        """
        Parameters:
            iterable (Sequence):    An iterable object which length can be accessed with len()
                                    (not an iterator).

            desc (str):             Description to put before the counter that describes the task.
                                    Defaults to "Iterating".

            end_desc (str):         Message to print when the task is finished. Defaults to "Done.".
        """
        self.iterable = iterable
        self.desc = desc if desc is not None else "Iterating"
        self.end_desc = end_desc if end_desc is not None else "Done."
    
    def __len__(self: tpe.Self) -> int:
        return len(self.iterable)
    
    def __iter__(self: tpe.Self) -> tp.Any:
        for idx, obj in enumerate(self.iterable, start=1):
            print(f"{self.desc}: {idx}/{len(self.iterable)}", end="\r", flush=True)
            yield obj

        print(f"\n{self.end_desc}")