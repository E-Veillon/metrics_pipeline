#!/usr/bin/python
"""
Functions relative to Periodic Table's (PT) elements properties.
"""

import re
from itertools import filterfalse
import functools as ft
import typing as tp
from collections.abc import Callable

from pymatgen.core import SiteCollection, Structure, Composition, Element, Species, DummySpecies
from pymatgen.core.periodic_table import ElementType
from pymatgen.io.cif import CifBlock

from .common_asserts import check_type


FormulaLike = str | tp.Sequence[str | int | Element]

ALL_ELT_Z_TO_SYMBOL = dict(enumerate((str(elt) for elt in Element), start=1))
"""Dict of {Z: symbol} of all 118 elements of the periodic table."""

ALL_ELT_SYMBOL_TO_Z = dict([(symbol, z) for z, symbol in enumerate((str(elt) for elt in Element), start=1)])
"""Dict of {symbol: Z} of all 118 elements of the periodic table."""

ELEMENT_TUPLE: dict[str | int, tuple[str, int]] = (
    {z: (symbol, z) for z, symbol in enumerate((str(elt) for elt in Element), start=1)} |
    {symbol: (symbol, z) for z, symbol in enumerate((str(elt) for elt in Element), start=1)}
)
"""
Convenient dict to get the (symbol, Z) pair corresponding to any passed
valid chemical element symbol or Z.
"""

PMG_ELTS_CATEGORIES = {
    str(elt_grp.value) for elt_grp in ElementType if elt_grp != ElementType.quadrupolar
}
GRP_ELTS_CATEGORIES = {f"group_{i}" for i in range(1,18)}
PRD_ELTS_CATEGORIES = {f"period_{i}" for i in range(1,7)}
ALL_ELTS_CATEGORIES = PMG_ELTS_CATEGORIES | GRP_ELTS_CATEGORIES | PRD_ELTS_CATEGORIES


def get_elts_from_symbol_or_z(symbols_or_z: list[str | int]) -> dict[str, int]:
    """
    Get symbols and atomic numbers of specified elements symbol or atomic number.

    Parameters
    ----------
    symbols_or_z: list[str | int]
        List of element symbols or atomic numbers to parse.

    Returns
    -------
    dict[str, int]
        Dict of the form {"symbol": atomic number} containing specified elements.
        If an empty list is passed, an empty dict is returned.
    """
    if symbols_or_z == []:
        return {}
    
    elts_dict = {}
    for elt in symbols_or_z:
        elt = int(elt) if isinstance(elt, str) and elt.isdecimal() else elt
        try:
            symbol, z = ELEMENT_TUPLE[elt]
        except KeyError:
            raise ValueError(f"Given symbol or atomic number {elt!r} is not a valid element.")

        elts_dict[symbol] = z

    return elts_dict


def get_elts_in_categories(cdts_list: list[str]) -> dict[str, int]:
    """
    Get symbols and atomic numbers of elements that are part of given categories.

    Parameters
    ----------
    cdts_list: list[str]
        List of string categories to parse.

    Returns
    -------
    dict[str, int]
        Dict of the form {"symbol": atomic number} containing elements that are parts of at least
        one of given categories. If an empty list is passed, an empty dict is returned.
    """
    if cdts_list == []:
        return {}

    # Parse conditions
    pmg_cdts = list(filter(lambda cdt: cdt in PMG_ELTS_CATEGORIES, cdts_list))
    grps_set = {
        int(cdt.split("_")[1])
        for cdt in filter(lambda cdt: cdt in GRP_ELTS_CATEGORIES, cdts_list)
    }
    prds_set = {
        int(cdt.split("_")[1])
        for cdt in filter(lambda cdt: cdt in PRD_ELTS_CATEGORIES, cdts_list)
    }
    # Embed all conditions in a single lambda
    predicate_func: Callable[[Element], bool] = lambda elt: (
        any(getattr(elt, f"is_{cdt}")() for cdt in pmg_cdts) or
        elt.group in grps_set or elt.row in prds_set
    )
    # Iterate on all elements
    all_elts = (Element(symbol) for symbol in ALL_ELT_SYMBOL_TO_Z)
    found_elts: dict[str, int] = dict(
        [(elt.symbol, elt.Z) for elt in filter(predicate_func, all_elts)]
    )
    return found_elts


def has_elements(
    structure: Structure | str, elements: list[str], format: str | None = None
) -> bool:
    """
    Whether the structure data contains at least one of given elements.

    Parameters
    ----------
    struccture: Structure | str
        The structure data to search into.

    elements: list[str]
        List of element symbols to search for.

    format: str, optional
        If `data` is passed as a str, precise the formatting ("cif" or "poscar").

    Returns
    -------
    bool
        `True` if the data contains at least one of the elements, else `False`.
    """
    if isinstance(structure, Structure):
        return any(elt in elements for elt in structure.composition.get_el_amt_dict().keys())

    if not isinstance(structure, str):
        raise TypeError(
            f"'data' expected a type 'Structure' or 'str', got {type(structure).__name__!r}."
        )
    if not elements:
        return False


    match format:
        case "cif":
            # Search for any label beginning with '_chemical_formula' containing composition
            formula_line = re.compile(fr"^_chemical_formula[a-zA-Z0-9_]+\s+(.+)$", flags=re.MULTILINE)
            for match in re.finditer(formula_line, structure):
                # Extract formula from the line
                formula = match.group(1)
                # Next line if no formula in the line
                if not formula:
                    continue
                # Parse element symbols from formula
                elts = re.findall(r"[A-Z][a-z]?", formula)
                if any(elt in elements for elt in elts):
                    return True
            else:
                # None of the formula lines contained any of searched elements
                return False
        case "poscar":
            # Search in the POSCAR element line (VASP 5.0+ format)
            elts = structure.split("\n", maxsplit=6)[5].split()
            return any(elt in elements for elt in elts)
        case str():
            raise ValueError(f"Unsupported format {format!r}.")
        case None:
            raise ValueError("'format' must be given when passing data as a string.")
        case _:
            raise TypeError(f"'format' expected a type 'str', got {type(format).__name__!r}.")


@tp.overload
def filter_by_elements(
    structures: list[Structure], elements: list[str], format: str | None = None
) -> tuple[list[Structure], int]: ...

@tp.overload
def filter_by_elements(
    structures: list[str], elements: list[str], format: str | None = None
) -> tuple[list[str], int]: ...

def filter_by_elements(
    structures: list[Structure] | list[str], elements: list[str], format: str | None = None
) -> tuple[list[Structure] | list[str], int]:
    """
    Filter structures containing at least one of given elements.

    Parameters
    ----------
    structures: list[Structure] | list[str]
        List of structure data to filter.

    elements: list[str]
        List of element symbols to search for.

    format: str, optional
        If structure data are passed as strings, precise the formatting ("cif" or "poscar").

    Returns
    -------
    list[Structure] | list[str]
        The filtered list of structure data that do not contain any searched element.
    int
        The number of filtered data.
    """
    has_elts = ft.partial(has_elements, elements=elements, format=format)
    filtered_structs = list(filterfalse(has_elts, structures))
    nbr_discarded = len(structures) - len(filtered_structs)
    return filtered_structs, nbr_discarded


def has_rare_gas(structure: SiteCollection | str) -> bool:
    """
    Searches for rare gas symbols in structural formula.

    Parameters:
        structure (SiteCollection | str): A pymatgen structure or a chemical formula.

    Returns:
        bool: True if the formula contains rare gases, False otherwise.
    """
    check_type(structure, "structure", (SiteCollection, str))

    if isinstance(structure, SiteCollection):
        return structure.composition.contains_element_type("noble_gas")

    return re.search(r"(He|Ne|Ar|Kr|Xe|Rn|Og)", structure) is not None


def discard_rare_gas_structures(
        structures: tp.Sequence[SiteCollection | str]
    ) -> tuple[list[SiteCollection | str], int]:
    """
    Eliminates structures containing rare gases and counts the number eliminated.

    Parameters:
        structures ([SiteCollection | str]): the structure data to scan.
    
    Returns:
        List[SiteCollection | str]: The list of data not containing rare gases.
        Int: The number of structures discarded.
    """
    check_type(structures, "structures", (tp.Sequence,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (SiteCollection, str))

    structures = list(structures)
    if not structures:
        return [], 0

    nbr_discarded = len(list(filter(has_rare_gas, structures)))
    kept_structs  = list(filterfalse(has_rare_gas, structures))

    return kept_structs, nbr_discarded


def has_rare_earth(structure: SiteCollection | str) -> bool:
    """
    Searches for rare earth (f block elements) symbols in structural formula.

    Parameters:
        structure (SiteCollection | str): A pymatgen structure or a chemical formula.

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


def discard_rare_earth_structures(
        structures: tp.Sequence[SiteCollection | str]
    ) -> tuple[list[SiteCollection | str], int]:
    """
    Eliminates structures containing rare earth elements and counts the number eliminated.

    Parameters:
        structures ([SiteCollection | str]): the structure data to scan.
    
    Returns:
        List[SiteCollection | str]: The list of data not containing rare earth elements.
        Int: The number of structures discarded.
    """
    check_type(structures, "structures", (tp.Sequence,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (SiteCollection, str))

    structures = list(structures)
    if not structures:
        return [], 0

    nbr_discarded = len(list(filter(has_rare_earth, structures)))
    kept_structs  = list(filterfalse(has_rare_earth, structures))

    return kept_structs, nbr_discarded


def get_elements(elts_data: FormulaLike) -> list[Element|Species|DummySpecies]:
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
    check_type(elts_data, "elts_data", (str, tp.Sequence))

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


def get_elemental_subsets(
        main_elts_set: FormulaLike,
        elts_subsets: tp.Sequence[FormulaLike]
    ) -> list[FormulaLike]:
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
    check_type(elts_subsets, "elts_subsets", (str, tp.Sequence))
    for idx, data in enumerate(elts_subsets):
        check_type(data, f"elts_data[{idx}]", (str, int, Element))

    ref_elts = get_elements(main_elts_set)

    sub_pd_list = list(filter(
        lambda pd_elts: all(elt in ref_elts for elt in get_elements(pd_elts)),
        elts_subsets
    ))

    return sub_pd_list


def get_element_group(
        atom: Element | str,
        return_type: tp.Literal["int", "str"] = "str"
    ) -> int | str:
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


def get_all_elements_groups(structure: SiteCollection | str) -> list[str]:
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


def get_element_valence_electrons(atom: str | Element) -> int:
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
