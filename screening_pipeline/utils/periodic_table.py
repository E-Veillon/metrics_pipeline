'''
Functions relative to Periodic Table's (PT) elements properties.
'''


########################################
# TYPE HINTING

from typing import Union, Iterable, List, Tuple, Literal, Sequence

########################################
# OPTIMIZATION MODULES

import sys
import re
from itertools import filterfalse

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import SiteCollection, Composition
from pymatgen.core.periodic_table import Element
from pymatgen.io.cif import CifBlock

########################################
# LOCAL MODULES

#from screening_pipeline.utils.fitted_values import EL_PER_XC_VOL
from screening_pipeline.utils.custom_types import FormulaLike

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

    assert isinstance(structures, Iterable), \
    f"Provided 'structures' argument is not iterable (got {type(structures)} instead)."

    structures = list(structures)
    if not len(structures): return [], 0

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

    assert isinstance(structures, Iterable), \
    f"Provided 'structures' argument is not iterable (got {type(structures)} instead)."

    structures = list(structures)
    if not len(structures): return [], 0

    nbr_discarded = len(list(filter(has_rare_earth, structures)))
    kept_structs  = list(filterfalse(has_rare_earth, structures))
    
    return kept_structs, nbr_discarded

def get_elements(
        elts_data: Union[str, Sequence[Union[str,int,Element]]]
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

    assert isinstance(elts_data, Iterable), \
    f"Provided 'elts_data' argument is not iterable (got {type(elts_data)} instead)."

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

def get_element_group(
        atom: Union[Element, str], 
        return_type: Literal['int', 'str'] = 'str'
    ) -> Union[int, str]:
    '''
    Finds the Periodic Table group of an element given as a string or Element object.

    Parameters:
        atom (str|Element):         The element to find the group for.

        return_type ('int'|'str'):  Whether to return the group number as an int
                                    (from 1 to 18, returns 19 for Lanthanides and
                                    20 for Actinides), or as a str representing the 
                                    group relative to electronic structure ('S1', 
                                    'S2', then from 'D1' to 'D10', then from 'P1' 
                                    to 'P6', 'L' for Lanthanides and 'A' for Actinides).
    
    Returns:
        int|str:    For return_type = 'str', the group of the element in 
                    '[block][group number in block]' format, e.g. for Fe it will return 'D6'.
                    For f-block elements, a lone letter will be returned, as Δ-Sol method
                    is not usable for these at the moment.
                    For return_type = 'int', the number of the group of the element, 
                    Lanthanides considered in "group 19" and Actinides in "group 20" as to
                    have a unique return for each group of elements, even if according to 
                    periodic table the best classification should be group 3 for both.
    '''

    assert return_type == 'str' or return_type == 'int'

    try: atom = Element(atom)
    except TypeError:
        raise TypeError(
            f"'atom' arg expected a 'str' or 'Element' type, got {type(atom)} instead."
        )
    except ValueError:
        raise ValueError(
            f"'atom' arg value '{atom}' is not recognized as an element."
        )

    if atom.block == 's':
        return (
            (atom.block.upper() + str(atom.group)) 
            if return_type == 'str' else atom.group
        )
    elif atom.block == 'p':
        return (
            (atom.block.upper() + str(atom.group - 12)) 
            if return_type == 'str' else atom.group
        )
    elif atom.block == 'd':
        return (
            (atom.block.upper() + str(atom.group - 2)) 
            if return_type == 'str' else atom.group
        )
    elif atom.block == 'f':
        return (
            ('L' if return_type == 'str' else 19) 
            if atom.row == 6 else ('A' if return_type == 'str' else 20)
        )
    else: raise NotImplementedError(f"{atom.block}: Unrecognized Periodic Table block.")

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
    
    elif isinstance(structure, str):
        formula   = CifBlock.from_str(structure).data["_chemical_formula_structural"]
        comp      = Composition(formula)
        elts_list = list(comp.keys())

    grps_list = list(map(get_element_group, elts_list))
    return grps_list

def get_element_valence_electrons(atom: Union[str,Element]) -> int:
    '''
    Gets the number of valence electrons of an element according to its group.

    Parameters:
        atom (str|Element): The element to compute the number of valence electrons from.
    
    Returns:
        int: The number of valence electrons corresponding to the element's group.
    '''
    group = get_element_group(atom, return_type='str')

    if group.startswith("S"):
        nb_val_elec = int(group[1])

    elif group.startswith("D") or group.startswith("P"):
        nb_val_elec = int(group[1]) + 2
        # Δ-Sol method counts all outermost s and d electrons in transition metals, 
        # even for d10 ones.
    elif group == 'L' or group == 'A':
        raise NotImplementedError(
            'f-block elements are not yet supported.'
        )
    else: raise ValueError('Provided string is not a recognized element group.')

    return nb_val_elec

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
        nbr_val_elec += get_element_valence_electrons(elt) * int(number)
    
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
        if elt.block == 's' or elt.block == 'p':
            continue
        elif elt.block == 'd': 
            val_elec_type = 'spd'
            break
        elif elt.block == 'f':
            raise NotImplementedError(
                'f-block elements are not supported in Δ-Sol method.'
            )
        else:
            raise ValueError(
                'Something is wrong with this function or Element objects "block" property.'
            )
    
    N_0        = get_all_valence_electrons(structure)
    value_name = '_'.join((dft_functional, val_elec_type))
    N_star     = EL_PER_XC_VOL[n_star_type][value_name]
    n          = float(N_0) / float(N_star)

    return n
