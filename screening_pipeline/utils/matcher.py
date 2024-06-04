"""
Functions to sort structures by stoichiometry and discard duplicates.
"""


########################################
# TYPE HINTING

from typing import Tuple, List, Union, Sequence

########################################
# OPTIMIZATION MODULES

import itertools
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure, SiteCollection
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import Entry
from pymatgen.analysis.structure_matcher import StructureMatcher

########################################
# LOCAL MODULES

from screening_pipeline.utils.utils import flatten

########################################
# LOCAL FUNCTIONS


def hash_stoichiometry(comp: Union[Composition, Entry, SiteCollection]) -> int:
    """
    Generate a hash from the fractional composition of a compatible pymatgen object, 
    ie. a Composition object or an object having a .composition attribute returning a
    Composition object.

    Parameters:
        comp (Entry|SiteCollection|Composition):    A pymatgen object containing a 
                                                    composition formula.

    Returns:
        int: The hash.
    """

    assert isinstance(comp, (Composition, Entry, SiteCollection)), \
    f"Given object type is not supported ({type(comp)})."
    
    if isinstance(comp, Composition):
        return hash(comp.fractional_composition)
    
    return hash(comp.composition.fractional_composition)

########################################

def group_by_stoichiometry(
        comps: Sequence[Union[Composition, Entry, SiteCollection]]
        ) -> List[List[Union[Composition, Entry, SiteCollection]]]:
    """
    Group Composition objects or objects having a .composition attribute by fractional 
    composition using the `hash_stoichiometry` hash function. Note that objects from 
    different classes but having their respective associated composition equal will be
    grouped together anyway.

    Parameters:
        comps ([Composition|Entry|SiteCollection]): The sequence of objects to group by
                                                    their composition.

    Returns:
        A list containing lists of objects with the same fractionnal composition.
    """

    assert isinstance(comps, Sequence)
    if not comps:
        return []

    sorted_comps = sorted(comps, key=hash_stoichiometry)

    return [
        list(grouped)
        for _, grouped in itertools.groupby(sorted_comps, hash_stoichiometry)
    ]

########################################

def hash_composition(comp: Union[Composition, Entry, SiteCollection]) -> int:
    """
    Generate a hash from the fractional composition of a compatible pymatgen object, 
    ie. a Composition object or an object having a .composition attribute returning a
    Composition object.

    Parameters:
        comp (Entry|SiteCollection|Composition):    A pymatgen object containing a 
                                                    composition formula.

    Returns:
        int: The hash.
    """

    assert isinstance(comp, (Composition, Entry, SiteCollection)), \
    f"Given object type is not supported ({type(comp)})."
    
    if isinstance(comp, Composition):
        return hash(comp)
    
    return hash(comp.composition)

########################################

def group_by_composition(
        comps: Sequence[Union[Composition, Entry, SiteCollection]]
        ) -> List[List[Union[Composition, Entry, SiteCollection]]]:
    """
    Group Composition objects or objects having a .composition attribute by
    their contained element types. Note that objects from different classes
    but having their respective associated composition equal will be grouped
    together anyway.

    Parameters:
        comps ([Composition|Entry|SiteCollection]): The sequence of objects to group by
                                                    their composition.

    Returns:
        A list of lists of objects containing the same elements.
    """

    assert isinstance(comps, Sequence)
    if not comps:
        return []

    sorted_comps = sorted(comps, key=hash_composition)

    return [
        list(grouped)
        for _, grouped in itertools.groupby(sorted_comps, hash_composition)
    ]

########################################

def _group_by_equivalence(structures: List[Structure]) -> List[List[Structure]]:
    """
    Group structures by equivalence using the StructureMatcher object.

    Parameters:
        structure (List[Structure]): The list of structure to match.

    Returns:
        List[List[Structure]]: A list containing lists of equivalent structures.
    """

    matcher = StructureMatcher()
    return matcher.group_structures(structures)

########################################

def batch_group_by_equivalence(
        structures: Sequence[Structure],
        workers: int = 1,
        comment: str = None
    ) -> List[List[List[Structure]]]:
    """
    Group structures by equivalence in two steps:
    First, groups by stoichiometry, then pass each sub-group in
    the pymatgen StructureMatcher in parallel for efficiency.


    Parameters:
        structures ([Structure]):   The list of structures to match.

        workers (int):              Number of parallel processes to spawn.
                                    Defaults to 1.

        comment (str):              Optional message to print next to tqdm
                                    progression bar.

    Returns:
        List[List[List[Structure]]]: A list containing lists of same composition
        containing lists of equivalent structures.
    """

    nbr_struct    = len(structures)
    chunksize     = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    grouped_structs = group_by_stoichiometry(structures)

    equivalent_structs = process_map(
        _group_by_equivalence,
        grouped_structs,
        max_workers=workers,
        chunksize=chunksize,
        desc=comment
    )

    return equivalent_structs


########################################


def remove_equivalent(
    structures: List[Structure], 
    workers: int = 1, 
    keep_equivalent: bool = False
) -> Tuple[List[Structure], int]:
    """
    Group structures by equivalence using multiple processes, then discards the duplicates.

    Parameters:
        structures (List[Structure]): The list of structures to match.

        workers (int):                The number of parallel processes to use.
                                      Defaults to 1.

        keep_equivalent (bool):       Whether to keep equivalent structures.
                                      If True, structures will be sorted by equivalence but
                                      not be discarded. Defaults to False.

    Returns:
        List[Structure]: The list of unique (or sorted) structures.
        Int: The number of discarded structures.
    """

    process_description = ("removing duplicates" if not keep_equivalent else "sorting structures")

    equivalent_structs = batch_group_by_equivalence(
        structures=structures, workers=workers, comment=process_description
    )

    nbr_discarded = 0

    if keep_equivalent:
        sorted_structs = flatten(equivalent_structs, level_of_flattening=2)
        return sorted_structs, nbr_discarded
    
    sorted_structs = flatten(equivalent_structs, level_of_flattening=1)
    nbr_discarded  = sum([len(sublist) - 1 for sublist in sorted_structs])
    unique_structs = [sublist[0] for sublist in sorted_structs]
    return unique_structs, nbr_discarded


########################################


def _get_novel_structures(
        structures: List[List[Structure]],
        dataset: List[Structure]
):
    """
    Filter out structures that are neither coming from the given dataset nor 
    equivalent to one of them.

    Parameters:
        structures ([[[Structure]]]):   Nested list of structures as the output from
                                        batch_group_by_equivalence().

        dataset ([Structure]):          Reference dataset of non-novel structures.
    
    Returns: List[Structure]
    The list of novel structures not seen in the dataset.
    """
    return flatten(
        list(
            filter(
                lambda l: len(l) == 1 and l[0] not in dataset,
                structures
            )
        )
    )
#----------------------------------------
def batch_get_novel_structures(
        structures: List[List[List[Structure]]],
        dataset: List[Structure],
        workers: int = 1
    ) -> List[Structure]:
    """
    Filter out structures that are neither coming from the given dataset nor 
    equivalent to one of them. Can be parallelized over compositional lists.

    Parameters:
        structures ([[[Structure]]]):   Nested list of structures as the output from
                                        batch_group_by_equivalence().

        dataset ([Structure]):          Reference dataset of non-novel structures.

        workers (int):                  Number of parallel processes to spawn.
    
    Returns: List[Structure]
    The list of novel structures not seen in the dataset.
    """
    if not isinstance(workers, int):
        raise TypeError(
            f"'workers arg expected a type 'int', got {type(workers)} instead."
        )
    if not workers > 0:
        raise ValueError(
            f"'workers' arg must be strictly positive (got {workers})."
        )
    
    get_novel_structs = partial(_get_novel_structures, dataset=dataset)

    novel_structs = process_map(
        get_novel_structs,
        structures,
        max_workers=workers,
        desc="search for novel structures"
    )
    novel_structs = flatten(novel_structs)
    return novel_structs
