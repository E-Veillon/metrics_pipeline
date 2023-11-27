# TYPE HINTING

from typing import Tuple, List, Union, Iterable

########################################
# OPTIMIZATION MODULES

import itertools
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS

from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

########################################


def hash_stoichiometry(structure: Structure) -> int:
    '''
    Generate a hash from the fractional composition of a structure.

    Parameters:
        structure (Structure): A pymatgen Structure object.

    Returns:
        int: The hash.
    '''

    return hash(structure.composition.fractional_composition)


def group_by_stoichiometry(structures: List[Structure]) -> List[List[Structure]]:
    '''
    Group structures by fractional composition using the `hash_stoichiometry` hash function.

    Parameters:
        structure (List[Structure]): The list of structures to group.

    Returns:
        List[List[Structure]]: A list containing lists of structures with the same fractional composition.
    '''

    sorted_structs = sorted(structures, key=hash_stoichiometry)
    return [
        list(grouped)
        for _, grouped in itertools.groupby(sorted_structs, hash_stoichiometry)
    ]


def group_by_equivalence(structures: List[Structure]) -> List[List[Structure]]:
    '''
    Group structure by equivalence using the StructureMatcher object.

    Parameters:
        structure (List[Structure]): The list of structure to match.

    Returns:
        List[List[Structure]]: A list containing lists of equivalent structures.
    '''

    matcher = StructureMatcher()
    return matcher.group_structures(structures)

def flatten(iterable: Iterable, level_of_flattening: int = 1) -> List:
    '''
    Unpacks a nested list or tuple without modifying elements order.

    Parameters:
        iterable (Iterable):        The iterable to unpack.

        level_of_flattening (Int):  The number of nested levels to unpack.
                                    Defaults to 1.
    
    Returns:
        A list flattened the specified number of times.
    '''

    for _ in range(1, level_of_flattening + 1):
        iterable = list(itertools.chain.from_iterable(iterable))
    return iterable

def remove_equivalent(
    structures: List[Structure], 
    workers: int = 1, 
    keep_equivalent: bool = False
) -> Tuple[List[Structure], int]:
    '''
    Group structures by equivalence using multiple processes.

    Parameters:
        structures (List[Structure]): The list of structures to match.

        workers (int):                The number of parallel processes to use.
                                      Defaults to 1.

        keep_equivalent (bool):       Whether to keep equivalent structures.
                                      If True, structures will be sorted by equivalence but not be discarded.
                                      Defaults to False.

    Returns:
        List[Structure]: The list of unique structures (or sorted structures if keep_equivalent = True).
        Int: The number of discarded structures.
    '''

    nbr_struct    = len(structures)
    chunksize     = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
    nbr_discarded = 0

    grouped_structs = group_by_stoichiometry(structures)

    equivalent_structs = process_map(
        group_by_equivalence,
        grouped_structs,
        max_workers=workers,
        chunksize=chunksize,
        desc='removing equivalents',
    )

    if keep_equivalent:
        sorted_structs = flatten(equivalent_structs, 2)
        return sorted_structs, nbr_discarded
    
    sorted_structs = flatten(equivalent_structs, 1)
    nbr_discarded  = sum([len(sublist) - 1 for sublist in sorted_structs])
    unique_structs = [sublist[0] for sublist in sorted_structs]
    return unique_structs, nbr_discarded
