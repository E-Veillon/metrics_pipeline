'''
Functions to sort structures by stoichiometry and discard duplicates.
'''


########################################
# TYPE HINTING

from typing import Tuple, List, Union, Sequence

########################################
# OPTIMIZATION MODULES

import itertools
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure, SiteCollection
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import Entry
from pymatgen.analysis.structure_matcher import StructureMatcher

########################################


def hash_stoichiometry(comp: Union[Composition, Entry, SiteCollection]) -> int:
    '''
    Generate a hash from the fractional composition of a compatible pymatgen object, 
    ie. a Composition object or an object having a .composition attribute returning a
    Composition object.

    Parameters:
        comp (Entry|SiteCollection|Composition):    A pymatgen object containing a 
                                                    composition formula.

    Returns:
        int: The hash.
    '''

    assert isinstance(comp, (Composition, Entry, SiteCollection)), \
    f'Given object type is not supported ({type(comp)}).'
    
    if isinstance(comp, Composition):
        return hash(comp.fractional_composition)
    
    return hash(comp.composition.fractional_composition)


def group_by_stoichiometry(
        comps: Sequence[Union[Composition, Entry, SiteCollection]]
        ) -> List[List[Union[Composition, Entry, SiteCollection]]]:
    '''
    Group Composition objects or objects having a .composition attribute by fractional 
    composition using the `hash_stoichiometry` hash function. Note that objects from 
    different classes but having their respective associated composition equal will be
    grouped together anyway.

    Parameters:
        comps ([Composition|Entry|SiteCollection]): The sequence of objects to group by
                                                    their composition.

    Returns:
        A list containing lists of objects with the same fractionnal composition.
    '''

    sorted_comps = sorted(comps, key=hash_stoichiometry)

    return [
        list(grouped)
        for _, grouped in itertools.groupby(sorted_comps, hash_stoichiometry)
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

def flatten(sequence: Sequence, level_of_flattening: int = 1) -> List:
    '''
    Unpacks a nested sequence without modifying elements order.

    Parameters:
        iterable (Iterable):        The iterable to unpack.

        level_of_flattening (Int):  The number of nested levels to unpack.
                                    Defaults to 1.
    
    Returns:
        A list flattened the specified number of times.
    '''

    for _ in range(1, level_of_flattening + 1):
        sequence = list(itertools.chain.from_iterable(sequence))
    return sequence

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

    grouped_structs     = group_by_stoichiometry(structures)
    process_description = ('removing duplicates' if not keep_equivalent else 'sorting structures')

    equivalent_structs = process_map(
        group_by_equivalence,
        grouped_structs,
        max_workers=workers,
        chunksize=chunksize,
        desc=process_description
    )

    if keep_equivalent:
        sorted_structs = flatten(equivalent_structs, level_of_flattening=2)
        return sorted_structs, nbr_discarded
    
    sorted_structs = flatten(equivalent_structs, level_of_flattening=1)
    nbr_discarded  = sum([len(sublist) - 1 for sublist in sorted_structs])
    unique_structs = [sublist[0] for sublist in sorted_structs]
    return unique_structs, nbr_discarded
