'''
Functions to find and discard structures containing rare gases.
'''


########################################
# TYPE HINTING

from typing import Union, Iterable, List, Tuple

########################################
# OPTIMIZATION MODULES

import re
from itertools import filterfalse

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import SiteCollection

########################################


def has_rare_gas(structure: Union[SiteCollection, str]) -> bool:
    """
    Searches for rare gas symbols in structural formula

    Parameters:
        structure (Union[SiteCollection, str]): A pymatgen structure or a chemical formula.

    Returns:
        bool: True if the formula contains rare gases, False otherwise.
    """

    assert isinstance(structure, (SiteCollection, str))

    if isinstance(structure, SiteCollection):
        return structure.composition.contains_element_type("noble_gas")

    return re.search(r"(He|Ne|Ar|Kr|Xe|Rn)", structure) is not None

def discard_rare_gas_structures(
        structures: Iterable[Union[SiteCollection, str]]
        ) -> Tuple[List[Union[SiteCollection, str]], int]:
    '''
    Eliminates structures containing rare gases and counts the number eliminated.

    Parameters:
        structures (Iterable[SiteCollection | str]]): the structure data to scan.
    
    Returns:
        List[Union[SiteCollection, str]]: The list of data not containing rare gases.
        Int: The number of structures discarded.
    '''

    nbr_discarded = len(list(filter(has_rare_gas, structures)))
    kept_structs  = list(filterfalse(has_rare_gas, structures))
    return kept_structs, nbr_discarded