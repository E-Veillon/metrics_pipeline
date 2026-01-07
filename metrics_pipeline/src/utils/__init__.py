"""Helper functionalities for other packages."""

from .common_asserts import check_type, check_num_value, raise_or_warn
from .visual_iterator import VisualIterator
from .flattener import flatten
from .redirect import redirect_c_stdout, redirect_c_stderr
from .periodic_table import (
    ALL_ELT_SYMBOL_TO_Z, ALL_ELT_Z_TO_SYMBOL,
    has_rare_gas, has_rare_earth, discard_rare_gas_structures, discard_rare_earth_structures,
    get_elements, get_elemental_subsets, get_element_group, get_all_elements_groups,
    get_element_valence_electrons, get_all_valence_electrons
)
from .spg_data import PG_TO_SYSTEM, SPG_NUM_TO_PG
from .parse_args import parse_input_args

__all__ = [
    "check_type", "check_num_value",
    "VisualIterator",
    "flatten",
    "redirect_c_stdout", "redirect_c_stderr",
    "ALL_ELT_SYMBOL_TO_Z", "ALL_ELT_Z_TO_SYMBOL",
    "has_rare_gas", "has_rare_earth", "discard_rare_gas_structures", "discard_rare_earth_structures",
    "get_elements", "get_elemental_subsets", "get_element_group", "get_all_elements_groups",
    "get_element_valence_electrons", "get_all_valence_electrons",
    "PG_TO_SYSTEM", "SPG_NUM_TO_PG",
    "parse_input_args"
]