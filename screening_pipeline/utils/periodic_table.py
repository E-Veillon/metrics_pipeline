'''
Functions relative to Periodic Table's (PT) elements properties.
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
from pymatgen.core.periodic_table import Element

########################################


def has_rare_gas(structure: Union[SiteCollection, str]) -> bool:
    """
    Searches for rare gas symbols in structural formula.

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

def get_element_group(atom: str|Element) -> str:
    '''
    Finds PT group of a given element as string or Element object.

    Parameters:
        atom (str|Element): The element to find group for.
    
    Returns:
        str: The group of the element in '[block][group number in block]' form,
             e.g. for Fe it will return 'D6'.
             For f-block elements, a lone letter will be returned, as Δ-Sol method
             is not usable for these at the moment.
    '''
    match atom:
        case Element():
            match atom.block:
                case 's': return atom.block.upper() + str(atom.group)
                case 'p': return atom.block.upper() + str(atom.group - 12)
                case 'd': return atom.block.upper() + str(atom.group - 2)
                case 'f': return ('L' if atom.row == 6 else 'A')
        case ('H'|'Li'|'Na'|'K'|'Rb'|'Cs'|'Fr'): return 'S1'
        case ('He'|'Be'|'Mg'|'Ca'|'Sr'|'Ba'|'Ra'): return 'S2'
        case ('Sc'|'Y'|'Lu'|'Lr'): return 'D1'
        case ('Ti'|'Zr'|'Hf'|'Rf'): return 'D2'
        case ('V'|'Nb'|'Ta'|'Db'): return 'D3'
        case ('Cr'|'Mo'|'W'|'Sg'): return 'D4'
        case ('Mn'|'Tc'|'Re'|'Bh'): return 'D5'
        case ('Fe'|'Ru'|'Os'|'Hs'): return 'D6'
        case ('Co'|'Rh'|'Ir'|'Mt'): return 'D7'
        case ('Ni'|'Pd'|'Pt'|'Ds'): return 'D8'
        case ('Cu'|'Ag'|'Au'|'Rg'): return 'D9'
        case ('Zn'|'Cd'|'Hg'|'Cn'): return 'D10'
        case ('B'|'Al'|'Ga'|'In'|'Tl'|'Nh'): return 'P1'
        case ('C'|'Si'|'Ge'|'Sn'|'Pb'|'Fl'): return 'P2'
        case ('N'|'P'|'As'|'Sb'|'Bi'|'Mc'): return 'P3'
        case ('O'|'S'|'Se'|'Te'|'Po'|'Lv'): return 'P4'
        case ('F'|'Cl'|'Br'|'I'|'At'|'Ts'): return 'P5'
        case ('Ne'|'Ar'|'Kr'|'Xe'|'Rn'|'Og'): return 'P6'
        case ('La'|'Ce'|'Pr'|'Nd'|'Pm'|'Sm'|'Eu'|'Gd'|'Tb'|'Dy'|'Ho'|'Er'|'Tm'|'Yb'):
            return 'L'
        case ('Ac'|'Th'|'Pa'|'U'|'Np'|'Pu'|'Am'|'Cm'|'Bk'|'Cf'|'Es'|'Fm'|'Md'|'No'):
            return 'A'
        case str(): return None
        case _: raise TypeError(f'expected a str, got {type(atom)}.')

def get_element_valence_electrons(atom: str|Element) -> int:
    '''
    Gets the number of valence electrons of an element according to its group.

    Parameters:
        atom (str|Element): The element to compute the number of valence electrons from.
    
    Returns:
        int: The number of valence electrons corresponding to the element's group.
    '''
    match get_element_group(atom):
        case ('S1'): return 1
        case ('S2'): return 2
        case ('D1'|'P1'): return 3
        case ('D2'|'P2'): return 4
        case ('D3'|'P3'): return 5
        case ('D4'|'P4'): return 6
        case ('D5'|'P5'): return 7
        case ('D6'|'P6'): return 8
        case ('D7'): return 9
        case ('D8'): return 10
        case ('D9'): return 11
        case ('D10'): return 12
        # Δ-Sol method counts all outermost s and d electrons in transition metals, 
        # even for d10 ones.
        case ('L'|'A'): raise NotImplementedError(
            'f-block elements are not taken into account yet.'
            )
        case None: raise ValueError('Provided string is not a recognized element')

def get_all_valence_electrons(structure: SiteCollection) -> int:
    '''
    Computes the number of valence electrons per unit cell.

    Parameters:
        structure (SiteCollection): The structure to compute.
    
    Returns:
        int: The number of valence electrons in the unit cell.
    '''
    
    nbr_val_elec = 0
    elts_dict    = structure.composition.element_composition.as_dict()

    for elt, number in elts_dict.items():
        nbr_val_elec += get_element_valence_electrons(elt)*int(number)
    
    return nbr_val_elec