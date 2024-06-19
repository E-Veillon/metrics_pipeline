"""
Functions to check structures validity, compare them, and discard duplicates.
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

def check_interatomic_distances(
        structures: Sequence[SiteCollection], 
        valid_tol: float = 0.5
    ) -> Tuple[List[SiteCollection], int]:
    '''
    Checks interatomic distances with respect to a tolerance in angstroms 
    for all given structures, then returns the valid ones in a list.

    Parameters:
        structures ([SiteCollection]):  The structures to check.

        valid_tol (float):              Tolerance below which a distance between 2 sites 
                                        is considered invalid and is eliminatory for the
                                        structure. Defaults to 0.5 angstroms.

    Returns:
        List[SiteCollection]: list of valid structures.
        int: Number of invalid structures discarded.
    '''
    def _is_valid(structure: SiteCollection):
        return structure.is_valid(tol=valid_tol)
    
    assert isinstance(structures, Sequence)
    assert all(isinstance(structure, SiteCollection) for structure in structures)
    assert isinstance(valid_tol, float)

    valid_structs = list(filter(_is_valid, structures))
    nbr_discarded = len(structures) - len(valid_structs)

    return valid_structs, nbr_discarded


########################################


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

def _group_by_equivalence(
        structures: List[Structure], test_volume: bool = False
    ) -> Tuple[List[List[Structure]], int]:
    """
    Group structures by equivalence using the StructureMatcher object.

    Parameters:
        structure (List[Structure]):    The list of structures to match.

        test_volume (bool):             Whether to test if some structures have an
                                        unphysical volume inferior than 1 Angström^3,
                                        and prevents them to pass in the matcher,
                                        assuming their unicity. Such structures may
                                        cause problems in the StructureMatcher, but are
                                        unlikely to happen in general, so it is recommended
                                        to only set it to True if problems arose
                                        in a first run due to that. Defaults to False.

    Returns: Tuple[List[List[Structure]], int]:
    A list containing lists of equivalent structures, and count of unmatched structures
    (equal to zero if test_volume is set to False).
    """
    unmatchables = []
    unmatch_count = 0

    if test_volume:
        unmatchables = list(filter(lambda t: t[1].volume < 1, enumerate(structures)))
        
        for idx, _ in reversed(unmatchables):
            structures.pop(idx)

        unmatchables = [[t[1]] for t in unmatchables] # Assumed to be uniques
        
        if unmatchables:
            unmatch_count += len(unmatchables)

    matcher = StructureMatcher()
    return (matcher.group_structures(structures) + unmatchables, unmatch_count)

########################################

def batch_group_by_equivalence(
        structures: Sequence[Structure],
        test_volume: bool = False,
        workers: int = 1,
        comment: str = None
    ) -> Tuple[List[List[List[Structure]]], int]:
    """
    Group structures by equivalence in two steps:
    First, groups by stoichiometry, then pass each sub-group in
    the pymatgen StructureMatcher in parallel for efficiency.


    Parameters:
        structures ([Structure]):   The list of structures to match.

        test_volume (bool):         Whether to test if some structures have an
                                    unphysical volume inferior than 1 Angström^3,
                                    and prevents them to pass in the matcher,
                                    assuming their unicity. Such structures may
                                    cause problems in the StructureMatcher, but
                                    are unlikely to happen in general, so it is
                                    only recommended to set it to True if problems
                                    arose in a first run due to that.
                                    Defaults to False.

        workers (int):              Number of parallel processes to spawn.
                                    Defaults to 1.

        comment (str):              Optional message to print next to tqdm
                                    progression bar.

    Returns: Tuple[List[List[List[Structure]]], int]: 
    A list containing lists of same composition containing lists of equivalent structures,
    and total count of unmatched structures (equal to zero if test_volume is set to False).
    """

    nbr_struct    = len(structures)
    chunksize     = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    grouped_structs = group_by_stoichiometry(structures)
    equiv_matcher = partial(_group_by_equivalence, test_volume=test_volume)

    match_results = process_map(
        equiv_matcher,
        grouped_structs,
        max_workers=workers,
        chunksize=chunksize,
        desc=comment
    )
    equivalent_structs = [t[0] for t in match_results] #type: List[List[List[Structure]]]
    total_unmatch_count = sum([t[1] for t in match_results])

    return (equivalent_structs, total_unmatch_count)


########################################


def remove_equivalent(
    structures: List[Structure],
    workers: int = 1,
    test_volume: bool = False,
    keep_equivalent: bool = False
) -> Tuple[List[Structure], int, int]:
    """
    Group structures by equivalence using multiple processes, then discards the duplicates.

    Parameters:
        structures (List[Structure]):   The list of structures to match.

        test_volume (bool):             Whether to test if some structures have an
                                        unphysical volume inferior than 1 Angström^3,
                                        and prevents them to pass in the matcher,
                                        assuming their unicity. Such structures may
                                        cause problems in the StructureMatcher, but
                                        are unlikely to happen in general, so it is
                                        only recommended to set it to True if problems
                                        arose in a first run due to that.
                                        Defaults to False.

        workers (int):                  The number of parallel processes to use.
                                        Defaults to 1.

        keep_equivalent (bool):         Whether to keep equivalent structures.
                                        If True, structures will be sorted by equivalence
                                        but not be discarded. Defaults to False.

    Returns:
        List[Structure]: The list of unique (or sorted) structures.
        Int: The number of discarded structures.
        Int: The number of unmatched structures (zero if test_volume = False).
    """

    process_description = ("removing duplicates" if not keep_equivalent else "sorting structures")

    equivalent_structs, nbr_unmatched = batch_group_by_equivalence(
        structures=structures,
        workers=workers,
        test_volume=test_volume,
        comment=process_description
    )

    nbr_discarded = 0

    if keep_equivalent:
        sorted_structs = flatten(equivalent_structs, level_of_flattening=2)
        return sorted_structs, nbr_discarded, nbr_unmatched
    
    sorted_structs = flatten(equivalent_structs, level_of_flattening=1)
    nbr_discarded  = sum([len(sublist) - 1 for sublist in sorted_structs])
    unique_structs = [sublist[0] for sublist in sorted_structs]
    return unique_structs, nbr_discarded, nbr_unmatched


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
    
    Returns: List[Structure]:
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
