#!/usr/bin/python
"""
Legacy module for matching structures similarity.
"""

import itertools as itt
import functools as ft
import typing as tp

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.analysis.structure_matcher import StructureMatcher

from metrics_pipeline.src.utils import check_type, check_num_value, flatten, VisualIterator
from metrics_pipeline.src.metrics.hash_matcher import group_compositions


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

    grouped_structs = group_compositions(list(structures), by="formula")
    equiv_matcher = ft.partial(_group_by_equivalence, test_volume=test_volume)

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


def remove_equivalent(
    structures: list[Structure] | list[SymmetrizedStructure],
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
