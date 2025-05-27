#!/usr/bin/python
"""
Functions relative to Periodic Table's (PT) elements properties.
"""


import re
from itertools import filterfalse
from typing import Union, List, Tuple, Literal, Sequence

# PYTHON MATERIAL GENOMICS
from pymatgen.core import SiteCollection, Composition, Element, Species, DummySpecies
from pymatgen.io.cif import CifBlock

# LOCAL IMPORTS
from .common_asserts import check_type
from .custom_types import FormulaLike


########################################


def has_rare_gas(structure: Union[SiteCollection, str]) -> bool:
    """
    Searches for rare gas symbols in structural formula.

    Parameters:
        structure (Union[SiteCollection, str]): A pymatgen structure or a chemical formula.

    Returns:
        bool: True if the formula contains rare gases, False otherwise.
    """
    check_type(structure, "structure", (SiteCollection, str))

    if isinstance(structure, SiteCollection):
        return structure.composition.contains_element_type("noble_gas")

    return re.search(r"(He|Ne|Ar|Kr|Xe|Rn|Og)", structure) is not None


########################################


def discard_rare_gas_structures(
        structures: Sequence[Union[SiteCollection, str]]
    ) -> Tuple[List[Union[SiteCollection, str]], int]:
    """
    Eliminates structures containing rare gases and counts the number eliminated.

    Parameters:
        structures ([SiteCollection | str]): the structure data to scan.
    
    Returns:
        List[Union[SiteCollection, str]]: The list of data not containing rare gases.
        Int: The number of structures discarded.
    """
    check_type(structures, "structures", (Sequence,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (SiteCollection, str))

    structures = list(structures)
    if not structures:
        return [], 0

    nbr_discarded = len(list(filter(has_rare_gas, structures)))
    kept_structs  = list(filterfalse(has_rare_gas, structures))

    return kept_structs, nbr_discarded


########################################


def has_rare_earth(structure: Union[SiteCollection, str]) -> bool:
    """
    Searches for rare earth (f block elements) symbols in structural formula.

    Parameters:
        structure (Union[SiteCollection, str]): A pymatgen structure or a chemical formula.

    Returns:
        bool: True if the formula contains rare earth, False otherwise.
    """
    check_type(structure, "structure", (SiteCollection, str))

    if isinstance(structure, SiteCollection):
        return structure.composition.contains_element_type("f-block")

    elts_grps_list = get_all_elements_groups(structure)
    has_lanthanoid = "L" in elts_grps_list
    has_actinoid   = "A" in elts_grps_list

    return has_lanthanoid or has_actinoid


########################################


def discard_rare_earth_structures(
        structures: Sequence[Union[SiteCollection, str]]
    ) -> Tuple[List[Union[SiteCollection, str]], int]:
    """
    Eliminates structures containing rare earth elements and counts the number eliminated.

    Parameters:
        structures ([SiteCollection | str]): the structure data to scan.
    
    Returns:
        List[Union[SiteCollection, str]]: The list of data not containing rare earth elements.
        Int: The number of structures discarded.
    """
    check_type(structures, "structures", (Sequence,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (SiteCollection, str))

    structures = list(structures)
    if not structures:
        return [], 0

    nbr_discarded = len(list(filter(has_rare_earth, structures)))
    kept_structs  = list(filterfalse(has_rare_earth, structures))

    return kept_structs, nbr_discarded


########################################


def get_elements(elts_data: FormulaLike) -> List[Element|Species|DummySpecies]:
    """
    Flexible converter to get a list of unique Element objects from a single string or any 
    iterable providing valid element symbols, atomic numbers, Element objects, or a mixture 
    of the three.

    Parameters:
        elts_data (str|[str|int|Element]):  The data to parse Elements objects from.
                                            If a single string is provided, it can either 
                                            be a raw formula (eg. "FePO4") or a composition 
                                            string containing element symbols separated by 
                                            "-" (eg. "Fe-P-O").
                                            If an iterable is given, it can contain valid 
                                            element symbols, atomic numbers and/or Element 
                                            objects.

    Raises: 
        ValueError if some of the given data does not represents valid elements.

    Returns: 
        A list of parsed Element objects.
    """
    check_type(elts_data, "elts_data", (str, Sequence))

    if isinstance(elts_data, str):
        return Composition(
            "".join(elts_data.split(sep="-")), strict=True
        ).element_composition.elements

    else:
        for idx, data in enumerate(elts_data):
            check_type(data, f"elts_data[{idx}]", (str, int, Element))

        return Composition(
            [(elt, 1) for elt in elts_data], strict=True
        ).element_composition.elements


########################################


def get_elemental_subsets(
        main_elts_set: FormulaLike,
        elts_subsets: Sequence[FormulaLike]
    ) -> List[FormulaLike]:
    """
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
    """
    check_type(elts_subsets, "elts_subsets", (str, Sequence))
    for idx, data in enumerate(elts_subsets):
        check_type(data, f"elts_data[{idx}]", (str, int, Element))

    ref_elts = get_elements(main_elts_set)

    sub_pd_list = list(filter(
        lambda pd_elts: all(elt in ref_elts for elt in get_elements(pd_elts)),
        elts_subsets
    ))

    return sub_pd_list


########################################


def get_element_group(
        atom: Union[Element, str],
        return_type: Literal["int", "str"] = "str"
    ) -> Union[int, str]:
    """
    Finds the Periodic Table group of an element given as a string or Element object.

    Parameters:
        atom (str|Element):         The element to find the group for.

        return_type ("int"|"str"):  Whether to return the group number as an int
                                    (from 1 to 18, returns 19 for Lanthanides and
                                    20 for Actinides), or as a str representing the 
                                    group relative to electronic structure ("S1", 
                                    "S2", then from "D1" to "D10", then from "P1" 
                                    to "P6", "L" for Lanthanides and "A" for Actinides).
    
    Returns:
        int|str:    For return_type = "str", the group of the element in 
                    "[block][group number in block]" format, e.g. for Fe it will return "D6".
                    For f-block elements, a lone letter will be returned, as Δ-Sol method
                    is not usable for these at the moment.
                    For return_type = "int", the number of the group of the element, 
                    Lanthanides considered in "group 19" and Actinides in "group 20" as to
                    have a unique return for each group of elements, even if according to 
                    periodic table the best classification should be group 3 for both.
    """
    check_type(atom, "atom", (str, Element))
    assert return_type in {"str", "int"}

    try:
        atom = Element(atom)
    except ValueError as exc:
        raise ValueError(
            f"'atom' arg value '{atom}' is not recognized as an element."
        ) from exc

    match atom.block:
        case "s":
            return (
                (atom.block.upper() + str(atom.group))
                if return_type == "str" else atom.group
            )
        case "p":
            return (
                (atom.block.upper() + str(atom.group - 12))
                if return_type == "str" else atom.group
            )
        case "d":
            return (
                (atom.block.upper() + str(atom.group - 2))
                if return_type == "str" else atom.group
            )
        case "f":
            return (
                ("L" if return_type == "str" else 19) if atom.row == 6
                else ("A" if return_type == "str" else 20)
            )
        case _:
            raise NotImplementedError(
                f"{atom.block}: Unrecognized Periodic Table block."
            )


########################################


def get_all_elements_groups(structure: Union[SiteCollection, str]) -> List[str]:
    """
    Finds PT group of each element contained in a structure.
    Provided structure can either be a pymatgen SiteCollection object,
    or a CIF formatted string containing a structural formula.

    Parameters:
        formula (SiteCollection|str): The structure to search element groups in.
    
    Returns:
        List[str]:  The group of each element in "[block][group number in block]" format,
                    e.g. for Fe it will return "D6", in a list.
                    For f-block elements, a lone letter will be returned, as Δ-Sol method
                    is not usable for these at the moment.
    """
    check_type(structure, "structure", (SiteCollection, str))

    if isinstance(structure, SiteCollection):
        elts_list = list(structure.composition.keys())

    elif isinstance(structure, str):
        formula   = CifBlock.from_str(structure).data["_chemical_formula_structural"]
        comp      = Composition(formula)
        elts_list = list(comp.keys())

    grps_list = list(map(get_element_group, elts_list))
    return grps_list # type: ignore


########################################


def get_element_valence_electrons(atom: Union[str,Element]) -> int:
    """
    Gets the number of valence electrons of an element according to its group.

    Parameters:
        atom (str|Element): The element to compute the number of valence electrons from.
    
    Returns:
        int: The number of valence electrons corresponding to the element's group.
    """
    check_type(atom, "atom", (str, Element))

    group = get_element_group(atom, return_type="str")

    if group.startswith("S"): # type: ignore
        nb_val_elec = int(group[1]) # type: ignore

    elif group.startswith("D") or group.startswith("P"): # type: ignore
        nb_val_elec = int(group[1]) + 2 # type: ignore
        # Δ-Sol method counts all outermost s and d electrons in transition metals,
        # even for d10 ones.
    elif group in {"L", "A"}:
        raise NotImplementedError(
            "f-block elements are not yet supported."
        )
    else: raise ValueError("Provided string is not a recognized element group.")

    return nb_val_elec


########################################


def get_all_valence_electrons(structure: SiteCollection) -> float:
    """
    Computes the number of valence electrons per unit cell.

    Parameters:
        structure (SiteCollection): The structure to compute.
    
    Returns:
        int: The number of valence electrons in the unit cell.
    """
    check_type(structure, "structure", (SiteCollection,))

    nbr_val_elec = 0
    elts_dict = structure.composition.element_composition.as_dict()

    for elt, number in elts_dict.items():
        nbr_val_elec += get_element_valence_electrons(elt) * int(number)

    nbr_val_elec -= structure.charge #e.g. a charge of +1 means there is 1 less electron

    return nbr_val_elec


########################################
