import re
from typing import Union

from pymatgen.core.structure import SiteCollection


def has_rare_gas(structure: Union[SiteCollection, str]) -> bool:
    """
    searches for rare gas symbols in structural formula

    Args:
        structure (Union[SiteCollection, str]): A pymatgen structure or a chemical formula
    Returns
        bool: True if the formula contains rare gases, False otherwise.
    """

    assert isinstance(structure, (SiteCollection, str))

    if isinstance(structure, SiteCollection):
        return structure.composition.contains_element_type("noble_gas")

    return re.search(r"(He|Ne|Ar|Kr|Xe|Rn)", structure) is not None
