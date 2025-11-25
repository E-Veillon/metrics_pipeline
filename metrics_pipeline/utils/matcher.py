#!/usr/bin/python
"""
Functions to check structures validity, compare them, and discard duplicates.
"""

import re
import itertools as itt
from functools import partial
import typing as tp
from pathlib import Path
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
import numpy as np

# PYTHON MATERIAL GENOMICS
from pymatgen.core import Structure, SiteCollection, Composition
from pymatgen.analysis.phase_diagram import Entry
from pymatgen.analysis.structure_matcher import StructureMatcher

# LOCAL IMPORTS
from .common_asserts import check_type, check_num_value
from .custom_types import PathLike
from .visual_iterator import VisualIterator
from .flattener import flatten
from .fitted_values import get_clementi_radius


########################################


def get_min_dist(radii_pair: tuple[float, float], valid_tol: float = 0.5) -> float:
            min_radii_dist = sum(radii_pair) * 0.9 # 10% uncertainty tolerance
            return max(valid_tol, round(min_radii_dist, 3)) # Conversion to Angströms


def check_interatomic_distances(
    structures: tp.Sequence[Structure],
    valid_tol: tp.Literal["radii"] | float = 0.5
) -> tuple[list[Structure], int]:
    '''
    Checks interatomic distances with respect to a tolerance in angstroms 
    for all given structures, then returns the valid ones in a list.

    Parameters:
        structures ([SiteCollection]):  The structures to check.

        valid_tol ("radii" | float):    Tolerance below which a distance between 2 sites 
                                        is considered invalid and is eliminatory for the
                                        structure. Defaults to 0.5 angstroms. If "radii"
                                        is passed, a more sophisticated scheme using pairwise
                                        atomic radii as calculated by Clementi et al. will
                                        be used to check valid distances.

    Returns:
        List[SiteCollection]: list of valid structures.
        int: Number of invalid structures discarded.
    '''
    check_type(structures, "structures", (tp.Sequence,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))
    check_type(valid_tol, "valid_tol", (float,))

    def _is_valid(structure: SiteCollection, valid_tol: tp.Literal["radii"] | float) -> bool:
        if valid_tol == "radii":
            radii = [get_clementi_radius(elt.symbol) for elt in struct.species]
            min_dists = np.array(list(map(get_min_dist, itt.combinations(radii, r=2))))
            true_dists = struct.distance_matrix[np.triu_indices(len(struct), 1)]
            return np.all(np.subtract(true_dists, min_dists) >= 0.0).item()
        return structure.is_valid(tol=valid_tol)

    is_valid = partial(_is_valid, valid_tol=valid_tol)
    valid_structs = list(
        filter(
            is_valid,
            VisualIterator(structures, desc="Checking interatomic distances")
        )
    )
    nbr_discarded = len(structures) - len(valid_structs)

    return valid_structs, nbr_discarded


########################################


def hash_stoichiometry(comp: tp.Union[Composition, Entry, SiteCollection]) -> int:
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
    check_type(comp, "comp", (Composition, Entry, SiteCollection))

    if isinstance(comp, Composition):
        return hash(comp.reduced_formula)

    return hash(comp.composition.reduced_formula)

########################################

def group_by_stoichiometry(
    comps: tp.Sequence[tp.Union[Composition, Entry, SiteCollection]]
) -> list[list[tp.Union[Composition, Entry, SiteCollection]]]:
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
    check_type(comps, "comps", (tp.Sequence,))
    if not comps:
        return []
    for idx, comp in enumerate(comps):
        check_type(comp, f"comps[{idx}]", (Composition, Entry, SiteCollection,))

    sorted_comps = sorted(comps, key=hash_stoichiometry)

    return [
        list(grouped)
        for _, grouped in itt.groupby(sorted_comps, hash_stoichiometry)
    ]

########################################

def hash_composition(comp: tp.Union[Composition, Entry, SiteCollection]) -> int:
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
    check_type(comp, "comp", (Composition, Entry, SiteCollection))

    if isinstance(comp, Composition):
        return hash(comp)

    return hash(comp.composition)


########################################


def group_by_composition(
    comps: tp.Sequence[tp.Union[Composition, Entry, SiteCollection]]
) -> list[list[tp.Union[Composition, Entry, SiteCollection]]]:
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
    check_type(comps, "comps", (tp.Sequence,))
    if not comps:
        return []
    for idx, comp in enumerate(comps):
        check_type(comp, f"comps[{idx}]", (Composition, Entry, SiteCollection,))

    sorted_comps = sorted(comps, key=hash_composition)

    return [
        list(grouped)
        for _, grouped in itt.groupby(sorted_comps, hash_composition)
    ]


########################################


def _group_by_equivalence(
    structures: list[Structure], test_volume: bool = False
) -> tuple[list[list[Structure]], int]:
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
    structures: tp.Sequence[Structure],
    test_volume: bool = False,
    workers: int | None = None,
    comment: str | None = None,
    sequential: bool = False
) -> tuple[list[list[list[Structure]]], int]:
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
                                    If not given, tqdm.contrib.concurrent.process_map
                                    default is used.

        comment (str):              Optional message to print next to tqdm
                                    progression bar.

        sequential (bool):          Whether to use sequential for-loop instead of multiprocessing
                                    scheme. If set to True, the 'workers' arg is ignored.
                                    Defaults to False.

    Returns: Tuple[List[List[List[Structure]]], int]: 
        A list containing lists of same composition containing lists of equivalent structures,
        and total count of unmatched structures (equal to zero if test_volume is set to False).
    """
    check_type(structures, "structures", (tp.Sequence,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))
    check_type(test_volume, "test_volume", (bool,))
    if workers is not None:
        check_type(workers, "workers", (int,))
        check_num_value(workers, "workers", ">", 0)
    if comment is not None:
        check_type(comment, "comment", (str,))

    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    grouped_structs = group_by_stoichiometry(structures)
    equiv_matcher = partial(_group_by_equivalence, test_volume=test_volume)

    if sequential:
        match_results = [
            equiv_matcher(grp) for grp in VisualIterator(grouped_structs, desc=comment)
        ]

    else:
        match_results = process_map(
            equiv_matcher,
            grouped_structs,
            max_workers=workers,
            chunksize=chunksize,
            desc=comment
        )

    equivalent_structs = [t[0] for t in match_results] #type: list[list[list[Structure]]]
    total_unmatch_count = sum(t[1] for t in match_results)

    return (equivalent_structs, total_unmatch_count)


########################################


def remove_equivalent(
    structures: list[Structure],
    workers: int|None = None,
    test_volume: bool = False,
    keep_equivalent: bool = False,
    sequential: bool = False
) -> tuple[list[Structure], int, int]:
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
                                        If not given, tqdm.contrib.concurrent.process_map
                                        default is used.

        keep_equivalent (bool):         Whether to keep equivalent structures.
                                        If True, structures will be sorted by equivalence
                                        but not be discarded. Defaults to False.

        sequential (bool):              Whether to use sequential for-loop instead of multiprocessing
                                        scheme. If set to True, the 'workers' arg is ignored.
                                        Defaults to False.

    Returns:
        List[Structure]: The list of unique (or sorted) structures.
        Int: The number of discarded structures.
        Int: The number of unmatched structures (zero if test_volume = False).
    """
    check_type(structures, "structures", (tp.Sequence,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))
    check_type(test_volume, "test_volume", (bool,))
    check_type(keep_equivalent, "keep_equivalent", (bool,))
    if workers is not None:
        check_type(workers, "workers", (int,))
        check_num_value(workers, "workers", ">", 0)

    process_description = (
        "removing duplicates" if not keep_equivalent
        else "sorting structures"
    )

    equivalent_structs, nbr_unmatched = batch_group_by_equivalence(
        structures=structures,
        workers=workers,
        test_volume=test_volume,
        comment=process_description,
        sequential=sequential
    )

    nbr_discarded = 0

    if keep_equivalent:
        sorted_structs = flatten(equivalent_structs, level_of_flattening=2)
        return sorted_structs, nbr_discarded, nbr_unmatched

    sorted_structs = flatten(equivalent_structs, level_of_flattening=1)
    nbr_discarded  = sum(len(sublist) - 1 for sublist in sorted_structs)
    unique_structs = [sublist[0] for sublist in sorted_structs]
    return unique_structs, nbr_discarded, nbr_unmatched


########################################


def _get_novel_structures(
    structures: list[list[Structure]],
    dataset: list[Structure]
):
    """
    Filter out structures that are neither coming from the given dataset nor 
    equivalent to one of them.

    Parameters:
        structures ([[Structure]]): List of lists of equivalent structures like
                                    the output from batch_group_by_equivalence().

        dataset ([Structure]):      Reference dataset of non-novel structures.
    
    Returns: List[Structure]:
        The list of novel structures not seen in the dataset.
    """
    # Novel structures did not match with any structure from the dataset,
    # therefore all sub-lists containing novel structures are the ones not
    # containing any structure from dataset, no matter how many structures
    # are in the sub-list.
    return flatten(
        list(
            filter(
                lambda grp: all(
                    struct != ref for struct, ref in itt.product(grp, dataset)
                ), structures
            )
        )
    )
#----------------------------------------
def batch_get_novel_structures(
    structures: list[list[list[Structure]]],
    dataset: list[Structure],
    workers: int | None = None
) -> list[Structure]:
    """
    Filter out structures that are neither coming from the given dataset nor 
    equivalent to one of them. Can be parallelized over compositional lists.

    Parameters:
        structures ([[[Structure]]]):   Nested list of structures as the output from
                                        batch_group_by_equivalence().

        dataset ([Structure]):          Reference dataset of non-novel structures.

        workers (int):                  Number of parallel processes to spawn.
                                        If not given, tqdm.contrib.concurrent.process_map
                                        default is used.
    
    Returns: List[Structure]
        The list of novel structures not seen in the dataset.
    """
    if workers is not None:
        check_type(workers, "workers", (int,))
        check_num_value(workers, "workers", ">", 0)

    get_novel_structs = partial(_get_novel_structures, dataset=dataset)
    chunksize = (min(len(structures) // 100, 10) if len(structures) >= 200 else 1)

    novel_structs = process_map(
        get_novel_structs,
        structures,
        max_workers=workers,
        chunksize=chunksize,
        desc="search for novel structures"
    )
    novel_structs = flatten(novel_structs)
    return novel_structs


########################################


def is_struct_dir(path: PathLike) -> bool:
    """
    Checks whether given path is a directory having
    naming convention of a structure directory.
    """
    #NOTE: Here 'path' must be a Path object, do not change for os.path.
    path = Path(path)
    return (
        path.is_dir()
        and re.match(r"\A[0-9]+_[A-Za-z0-9\(\)]+\Z", path.name) is not None
    )


########################################


def match_struct_dirs(
    path: PathLike,
    indices: list[int] | None = None,
    match_all: bool = True,
    unique: bool = True,
    no_return: bool = False
) -> list[str]:
    """
    Match structure directories at given path.
    
    Parameters:
        path (str|Path):    Where to search for structure directories.

        indices ([int]):    Indices of wanted structure directories to match.
                            If not given, matches all directories conforming to structure
                            directory format. In this default behavior, match_all, unique
                            and no_return are ignored.

        match_all (bool):   Whether each given index should match at least one directory.
                            If True, raise an error when an index does not match any directory.
                            Defaults to True.

        unique (bool):      Whether each given index should match at most one directory.
                            If True, raise an error when an index matches several directories.
                            Defaults to True.

        no_return (bool):   Whether to only match structure directories without
                            returning the matching paths list for memory saving.
                            If True, an empty list is always returned.
                            Defaults to False.

    Raises:
        ValueError if match_all is True and an index does not match any structure directory.
        ValueError if unique is True and an index does match with several structure directories.

    Returns:
        [str]: List of paths of matching structure directories.
    """
    subtree = list(Path(path).iterdir())

    if indices is None:
        return sorted(list(map(str, filter(is_struct_dir, subtree))))

    all_matching_dirs = []

    for idx in sorted(indices):
        matching_dirs = list(
            filter(
                lambda f: f.name.startswith(f"{idx}_") and is_struct_dir(f),
                subtree
            )
        )

        if match_all and not matching_dirs:
            raise ValueError(
                f"Index '{idx}' do not match with any structure "
                f"directory in {path}."
            )
        if unique and len(matching_dirs) > 1:
            raise ValueError(
                f"Index '{idx}' matches with more than one structure directory, "
                "which should not be possible as structure indexation is unique "
                "inside a same batch. Please make sure that there are no parasite "
                f"directories in {path}, i.e. non-structure directories starting with "
                f"'{idx}_' or structure directories moved from other batches."
            )
        if not no_return:
            all_matching_dirs += sorted(matching_dirs)

    return list(map(str, all_matching_dirs))


########################################
