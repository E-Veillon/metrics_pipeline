"""Functions dealing with the processing and checking of GenMat preprocessed structure names."""

import re
import typing as tp

if tp.TYPE_CHECKING:
    from pymatgen.core import Structure

from .visual_iterator import VisualIterator


_ELT_WITH_INDEX_REGEX = r"[A-Z][a-z]?\d*"
_GENMAT_NAME_REGEX = fr"^(\d+)_(?:{_ELT_WITH_INDEX_REGEX}|\((?:{_ELT_WITH_INDEX_REGEX})+\)\d*)+$"


def is_genmat_name(name: str) -> bool:
    """
    Whether given name follows GenMat preprocessed structure naming convention.
    """
    return re.fullmatch(_GENMAT_NAME_REGEX, name) is not None


def check_genmat_name(name: str) -> None:
    """
    Verify that given name corrresponds to preprocessed structure naming convention.
    Raises a `ValueError` if it does not.
    """
    if is_genmat_name(name):
        return
    raise ValueError(
        f"{name!r} does not follow GenMat preprocessed structure naming conventions. "
        "Please make sure your data was preprocessed with the 'preprocess.py' "
        "script before going further."
    )


def get_genmat_name_idx(name: str) -> int:
    """
    Extract the index part of a GenMat preprocessed structure name.
    """
    check_genmat_name(name)
    match = re.match(_GENMAT_NAME_REGEX, name)
    assert match is not None
    return int(match.group(1))


def generate_genmat_names(structures: list[Structure]) -> list[Structure]:
    """
    Generate GenMat indexed names based on the position of each structure in the list.
    The new name is stored inside the 'header' key of structures 'properties' attribute.
    If another name already exists in 'header', it is moved to the '_original_header' key.
    """
    iterator: VisualIterator[tuple[int, Structure]] = VisualIterator.from_big_iterator(
        enumerate(structures), n_elts=len(structures),
        desc="Generating GenMat names", unit="generated", percent=True
    )
    for idx, structure in iterator:
        if old_header:=structure.properties.get('header'):
            structure.properties["_original_header"] = old_header
        structure.properties["header"] = f"{idx}_{structure.reduced_formula}"

    return structures
