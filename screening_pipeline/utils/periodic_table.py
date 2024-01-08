'''
Functions relative to Periodic Table's (PT) elements properties.
'''


########################################
# TYPE HINTING

from typing import Union, Iterable, List, Tuple, Literal, Sequence

########################################
# OPTIMIZATION MODULES

import re
from itertools import filterfalse

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import SiteCollection, Composition
from pymatgen.core.periodic_table import Element
from pymatgen.io.cif import CifBlock

########################################
# LOCAL MODULES

#from screening_pipeline.utils import EL_PER_XC_VOL

########################################
# TYPE ALIASES

FormulaLike = Union[str, Iterable[Union[str, int, Element]]]

########################################
# LOCAL FUNCTIONS

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

    return re.search(r"(He|Ne|Ar|Kr|Xe|Rn|Og)", structure) is not None

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

def has_rare_earth(structure: Union[SiteCollection, str]) -> bool:
    """
    Searches for rare earth (f block elements) symbols in structural formula.

    Parameters:
        structure (Union[SiteCollection, str]): A pymatgen structure or a chemical formula.

    Returns:
        bool: True if the formula contains rare earth, False otherwise.
    """

    assert isinstance(structure, (SiteCollection, str))

    if isinstance(structure, SiteCollection):
        return structure.composition.contains_element_type("f-block")
    
    elts_grps_list = get_all_elements_groups(structure)
    has_lanthanoid = 'L' in elts_grps_list
    has_actinoid   = 'A' in elts_grps_list

    return has_lanthanoid or has_actinoid

def discard_rare_earth_structures(
        structures: Iterable[Union[SiteCollection, str]]
    ) -> Tuple[List[Union[SiteCollection, str]], int]:
    '''
    Eliminates structures containing rare earth elements and counts the number eliminated.

    Parameters:
        structures (Iterable[SiteCollection | str]]): the structure data to scan.
    
    Returns:
        List[Union[SiteCollection, str]]: The list of data not containing rare earth elements.
        Int: The number of structures discarded.
    '''

    nbr_discarded = len(list(filter(has_rare_earth, structures)))
    kept_structs  = list(filterfalse(has_rare_earth, structures))
    return kept_structs, nbr_discarded

def get_elements(
        elts_data: str|Iterable[str|int|Element]
    ) -> List[Element]:
    '''
    Flexible converter to get a list of unique Element objects from a single string or any 
    iterable providing valid element symbols, atomic numbers, Element objects, or a mixture 
    of the three.

    Parameters:
        elts_data (str|[str|int|Element]):  The data to parse Elements objects from.
                                            If a single string is provided, it can either 
                                            be a raw formula (eg. 'FePO4') or a composition 
                                            string containing element symbols separated by 
                                            '-' (eg. 'Fe-P-O').
                                            If an iterable is given, it can contain valid 
                                            element symbols, atomic numbers and/or Element 
                                            objects.

    Raises: 
        ValueError if some of the given data does not represents valid elements.

    Returns: 
        A list of parsed Element objects.
    '''

    assert isinstance(elts_data, Iterable)

    if isinstance(elts_data, str):
        elts_list = Composition(''.join(elts_data.split(sep='-')), strict=True).elements
    
    else:
        assert all([isinstance(elt, (str, int, Element)) for elt in elts_data])
        elts_list = Composition([(elt, 1) for elt in elts_data], strict=True).elements

    return elts_list

def get_elemental_subsets(
        main_elts_set: FormulaLike, 
        elts_subsets: Sequence[FormulaLike]
    ) -> List[str]:
    '''
    Flexible function to extract all formulas from a given sequence that are fully made 
    of same elements as the given main formula. The atomic fractions are not taken into 
    account, only presence and absence of the elements are checked.

    Parameters:
        main_elts_set (str|Iterable):   The reference formula to get sub-formulas from.
                                        Only formulas containing only elements that are 
                                        present in this one will be returned.
        
        elts_subsets ([str|Iterable]):  The pool of formulas from which subformulas must
                                        be extracted.
                
    Returns:
        List of the formulas fully included in the main one.
    '''
    
    ref_elts = get_elements(main_elts_set)

    sub_pd_list = list(filter(
        lambda pd_elts: all([elt in ref_elts for elt in get_elements(pd_elts)]), 
        elts_subsets
    ))

    return sub_pd_list

def get_element_group(atom: Union[Element, str]) -> str:
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
        case _: raise TypeError(f'expected a Element or str, got {type(atom)}.')

def get_all_elements_groups(structure: Union[SiteCollection, str]) -> List[str]:
    '''
    Finds PT group of each element contained in a structure.
    Provided structure can either be a pymatgen SiteCollection object,
    or a CIF formatted string containing a structural formula.

    Parameters:
        formula (SiteCollection|str): The structure to search element groups in.
    
    Returns:
        List[str]:  The group of each element in '[block][group number in block]' format,
                    e.g. for Fe it will return 'D6', in a list.
                    For f-block elements, a lone letter will be returned, as Δ-Sol method
                    is not usable for these at the moment.
    '''

    assert isinstance(structure, (SiteCollection, str)), \
    'Provided structure must be a SiteCollection object or a CIF string'

    if isinstance(structure, SiteCollection):
        elts_list = list(structure.composition.keys())
    
    if isinstance(structure, str):
        formula   = CifBlock.from_str(structure).data["_chemical_formula_structural"]
        comp      = Composition(formula)
        elts_list = list(comp.keys())

    grps_list = list(map(get_element_group, elts_list))
    return grps_list

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
    
    nbr_val_elec -= structure.charge #e.g. a charge of +1 means there is 1 less electron
    
    return nbr_val_elec

def get_delta_sol_el_ratio(
        structure: SiteCollection, 
        dft_functional: Literal['LDA','PBE','AM05'] = 'PBE', 
        n_star_type: Literal['MIN', 'BEST', 'MAX'] = 'BEST'
    ) -> float:
    '''
    Computes n = N0/N* the electron ratio to add or remove from 
    the structure in the Δ-Sol method developped by Chan et al.

    Reference:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (reference 32 in screening_pipeline/Bibliography)
    '''
    from screening_pipeline.utils import EL_PER_XC_VOL
    
    val_elec_type = 'sp'

    for elt in structure.elements:
        match elt.block:
            case ('s'|'p'): continue
            case 'd': 
                val_elec_type = 'spd'
                break
            case 'f': raise NotImplementedError('f-block elements are not taken into accoount in Δ-Sol method.')
            case _: raise ValueError('Something is wrong with this loop or Element objects "block" property.')
    
    N_0        = get_all_valence_electrons(structure)
    value_name = '_'.join(dft_functional, val_elec_type)
    N_star     = EL_PER_XC_VOL[n_star_type][value_name]
    n          = float(N_0) / float(N_star)

    return n
