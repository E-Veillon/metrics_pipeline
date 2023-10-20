from typing import Tuple, List, Union
import itertools

from tqdm.contrib.concurrent import process_map

from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher


def hash_stoichiometry(structure: Structure) -> int:
    """
    Generate a hash from the fractional composition of a structure.

    Args:
        structure (Structure): A pymatgen structure.
    Returns
        int: The hash.
    """

    return hash(structure.composition.fractional_composition)


def group_by_stoichiometry(structs: List[Structure]) -> List[List[Structure]]:
    """
    Group structure by fractional composition using the `hash_stoichiometry` hash function.

    Args:
        structure (List[Structure]): The list of structure.
    Returns
        List[List[Structure]]: A list containing lists of structures with the same fractional composition.
    """

    sorted_structs = sorted(structs, key=hash_stoichiometry)
    return [
        list(grouped)
        for _, grouped in itertools.groupby(sorted_structs, hash_stoichiometry)
    ]


def group_by_equivalance(structs: List[Structure]) -> List[List[Structure]]:
    """
    Group structure by equivalance using the StructureMatcher object.

    Args:
        structure (List[Structure]): The list of structure.
    Returns
        List[List[Structure]]: A list containing lists of equivalent structures.
    """

    matcher = StructureMatcher()
    return matcher.group_structures(structs)


def remove_equivalent(
    structures: List[Structure], workers: int = 1, keep_equivalent: bool = False
) -> Union[List[Structure], Tuple[List[Structure], List[Structure]]]:
    """
    Group structures by equivalence using multiple workers.

    Args:
        structure (List[Structure]): The list of structure.
        workers (int): The number of worker used.
        keep_equivalent (bool): A flag to keep all equivalent structures.
    Returns
        List[Structure]: The list of unique structures.
    """

    if len(structures) >= 200:
        chunksize = min(len(structures) // 100, 10)
    else:
        chunksize = 1

    grouped_structs = group_by_stoichiometry(structures)

    equivalent_struct = process_map(
        group_by_equivalance,
        grouped_structs,
        max_workers=workers,
        chunksize=chunksize,
        desc="removing equivalent",
    )
    equivalent_struct = sum(equivalent_struct, [])

    kept_structs = [lst_structs[0] for lst_structs in equivalent_struct]

    if keep_equivalent:
        return kept_structs

    duplicated_struct = sum([lst_structs[1:] for lst_structs in equivalent_struct], [])

    return kept_structs, duplicated_struct
