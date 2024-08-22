#!/usr/bin/python
"""
Functions to initialize phase diagrams with energy convex hull,
process phase diagram entries, and compute energy above hull of entries.
"""


from typing import Dict, Union, Sequence, Any, List, Optional, Set, Tuple
from functools import partial
from tqdm.contrib.concurrent import process_map

# PYTHON MATERIAL GENOMICS
from pymatgen.core import Element, Composition
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

# LOCAL IMPORTS
from .common_asserts import check_type, check_num_value
from .custom_types import FormulaLike
from .flattener import flatten
from .matcher import group_by_composition
from .periodic_table import get_elements


########################################


def get_max_dim(entries: Sequence[PDEntry]) -> int:
    """Search for the maximum number of distinct elements in given entries."""
    check_type(entries, "entries", (Sequence,))
    for idx, entry in enumerate(entries):
        check_type(entry, f"entries[{idx}]", (PDEntry,))

    return max(len(entry.elements) for entry in entries)


########################################


def init_entries_from_dict(
    entries_dict: Dict[str, Dict],
    attribute: Optional[str] = None
) -> List[PDEntry]:
    """
    Convert structures data dict into a list of PDEntry objects
    that are compatible with phase diagrams.

    Parameters:
        entries_dict (dict):    Dict of dicts containing following structure data keys:
                                - entry_id (either dataset ID or generated name),
                                - composition (Composition object),
                                - final_energy (float).

        attribute (str):        Label to put as attribute in the entries.
    
    Returns:
        List[PDEntry]: The list of phase diagram entries.
    """
    check_type(entries_dict, "entries_dict", (Dict,))
    for idx, key in enumerate(entries_dict.keys()):
        check_type(key, f"entries_dict.keys()[{idx}]", (str,))
    for idx, val in enumerate(entries_dict.values()):
        check_type(val, f"entries_dict.values()[{idx}]", (Dict,))
    if attribute is not None:
        check_type(attribute, "attribute", (str,))

    print(f"Convert {len(entries_dict)} '{attribute}' structures to entries...")
    entries_list = [
        PDEntry(
            composition=entry["composition"],
            energy=entry["final_energy"],
            name=entry["entry_id"],
            attribute=attribute
        ) for entry in entries_dict.values()
    ]
    print("Conversion finished.")
    return entries_list


########################################


def filter_database_entries(
    entries: Dict[str, Any]|List[PDEntry],
    ref_elts: Optional[List[Element]] = None,
    max_dim: Optional[int] = None
) -> List[PDEntry]:
    """
    Filter out entries in database file that are useless for generated data.

    Parameters:
        entries (dict):         Entries to filter.

        ref_elts ([Element]):   List of useful Elements to keep in the dataset.
                                All structures containing other elements are discarded.

        max_dim (int):          Max dimension of entries to keep.

    Returns:
        List[PDEntry]: List of useful entries.
    """

    check_type(entries, "entries", (Dict, List))
    if max_dim is not None:
        check_type(max_dim, "max_dim", (int,))
    if ref_elts is not None:
        check_type(ref_elts, "ref_elts", (List,))

    print("Filtering reference dataset...")

    if isinstance(entries, Dict):
        entries = init_entries_from_dict(entries, attribute="ref_struct")

    if ref_elts is not None:
        entries = _get_relevant_entries(entries, ref_elts)

    if max_dim is not None:
        entries = list(
            filter(
                lambda entry: len(entry.elements) <= max_dim,
                entries
            )
        )
    print(f"Filtering done, {len(entries)} entries left.")
    return entries


########################################


def get_elements_from_entries(entries: Sequence[PDEntry]) -> List[Element]:
    """
    Get a list of all unique elements present in a sequence of PDEntry objects.

    Parameters:
        entries ([PDEntry]): The entries to get unique elements from.

    Returns:
        The list of found Element objects.
    """
    check_type(entries, "entries", (Sequence,))
    if not entries:
        return []
    for idx, entry in enumerate(entries):
        check_type(entry, f"entries[{idx}]", (PDEntry,))

    if len(entries) == 1:
        return entries[0].elements

    return list(set(flatten([entry.elements for entry in entries])))


########################################


def group_by_dim_and_comp(
    entries: Sequence[PDEntry],
    max_dim: Optional[int] = None,
    workers: int = 1
) -> List[List[List[PDEntry]]]:
    """
    Group entries in sublists according to the number of elements,
    then group entries in each sublist in subsublists according to
    their elemental composition. Can be parallelized over each dim.

    Parameters:
        entries ([PDEntry]):    entries to sort out.

        max_dim (int):          the highest number of elements an entry can have.
                                If some entries have higher dim than given max_dim,
                                they will be discarded. If max_dim is higher than the
                                actual highest dim in entries, empty sublists will be
                                put in the additional list indexes. If max_dim is not
                                given, it will be infered from given entries.

        workers (int):          Number of parallel processes to spawn.

    Returns:
        List[List[List[PDEntry]]]: List of lists of entries of same dimension sorted
        in sublists according to their composition.
    """
    check_type(entries, "entries", (Sequence,))
    if not entries:
        return []
    for idx, entry in enumerate(entries):
        check_type(entry, f"entries[{idx}]", (PDEntry,))
    if max_dim is not None:
        check_type(max_dim, "max_dim", (int,))
    check_type(workers, "workers", (int,))

    def _group_one_dim(entries: Sequence[PDEntry], dim: int):
        """Group all entries of a specific dimension by composition."""
        dim_group = list(filter(lambda entry: len(entry.elements) == dim, entries))
        grouped_entries = group_by_composition(dim_group)
        return grouped_entries

    initial_group = [[]]  # fill the index 0 to match indexes and entries dimension

    if max_dim is None:
        max_dim = get_max_dim(entries)

    grouper = partial(_group_one_dim, entries=entries)

    grouped_entries = process_map(
        grouper,
        list(range(1, max_dim + 1)),
        max_workers=workers,
        chunksize=1,
        desc="Grouping entries by dim and comp"
    )

    grouped_entries = initial_group + grouped_entries
    return grouped_entries


########################################


def get_sub_entries(
    main_entry: PDEntry, entry_pool: Sequence[PDEntry]
) -> List[PDEntry]:
    """
    Extract all PDEntry objects whose elemental composition are subsets of the
    composition of the main PDEntry object from a sequence of entries.

    Parameters:
        main_entry (PDEntry): The reference entry to search sub-entries from.

        entry_pool ([PDEntry]): The sequence of entries to search sub-entries in.

    Returns:
        The list of all found sub-entries relative to the main entry.
    """

    check_type(main_entry, "main_entry", (PDEntry,))
    check_type(entry_pool, "entry_pool", (Sequence,))
    if not entry_pool:
        return []
    for idx, entry in enumerate(entry_pool):
        check_type(entry, f"entry_pool[{idx}]", (PDEntry,))

    sub_entries = list(
        filter(
            lambda entry: all(elt in main_entry.elements for elt in entry.elements),
            entry_pool,
        )
    )
    return sub_entries


########################################


def _get_relevant_entries(
    entries: Sequence[PDEntry],
    ref_elts: Union[Sequence[Element], set[Element]]
) -> List[PDEntry]:
    """
    Extract entries that are only composed of given reference Elements.

    Parameters:
        entries ([PDEntry]):    Entries to filter out.

        ref_elts ([Elements]):  Reference Elements defining the wanted entry space.

    Returns:
        A list of filtered entries.
    """

    ref_entry = PDEntry(
        composition=Composition([(elt, 1) for elt in ref_elts]),
        energy=0.0,
        name="ref_entry",
    )

    return get_sub_entries(main_entry=ref_entry, entry_pool=entries)


########################################


def get_lacking_elts_entries(
    entries: Sequence[PDEntry],
    ref_elts: Union[Sequence[Element], set[Element]]
) -> List[PDEntry]:
    """
    Check elemental entries with respect to given reference elements,
    then initialize lacking elemental entries with an energy of 0.0 eV.

    Parameter:
        entries ([PDEntry]):    Entries to check elemental entries in.

        ref_elts ([Elements]):  Reference Elements defining the wanted entry space.

    Returns:
        A list of auto-defined elemental entries.
    """
    check_type(entries, "entries", (Sequence,))
    check_type(ref_elts, "ref_elts", (Sequence, Set))

    if not entries:
        user_elt_entries = []
    else:
        for idx, entry in enumerate(entries):
            check_type(entry, f"entries[{idx}]", (PDEntry,))

        user_elt_entries = list(
            filter(
                lambda entry: entry.is_element and entry.elements[0] in ref_elts,
                entries
            )
        )
    user_defined_elts = get_elements_from_entries(user_elt_entries)
    lacking_elts = list(
        filter(
            lambda elt: elt not in user_defined_elts,
            ref_elts
        )
    )
    if not lacking_elts: # No lacking entry
        return []

    auto_defined_elts_entries = [
        PDEntry(
            composition=Composition(str(elt)),
            energy=0.0,
            name=elt.symbol,
            attribute="element_ref",
        ) for elt in lacking_elts
    ]

    return auto_defined_elts_entries


########################################


def _process_pd_data(
    entries: Optional[Sequence[PDEntry]] = None,
    ref_elts: Optional[FormulaLike] = None
) -> Tuple[List[PDEntry], List[Element]]:
    """
    Process entries and reference elements necessary to build a phase diagram.

    Parameters:
        entries ([PDEntry]):                The entries that will be put into the PhaseDiagram.
                                            If there is no unary entries for some given ref_elts,
                                            default unary entries with energy = 0.0 eV are created.
                                            Entries containing elements that are not referenced in
                                            ref_elts are ignored. This behaviour is particularly
                                            useful if one needs to initialize several diagrams
                                            from different parts of the same entry dataset.

        ref_elts (str|[str|int|Element]):   The reference elements of the phase diagram.
                                            Can be given as a single formula string
                                            (e.g. 'FePO4') or a composition string (eg. 'Fe-P-O'),
                                            or as a sequence containing valid element symbols,
                                            atomic numbers and/or Element objects.
                                            If not given, they are inferred from given entries.
                                            In that case, returned ref_elts contains needed
                                            Element objects to make the minimal phase diagram
                                            containing all given entries.

    Returns:
        Processed entries and elements to build a PhaseDiagram.
    """
    match (bool(entries), bool(ref_elts)):
        case (False, False): # Both 'entries' and 'ref_elts' are empty or None
            raise ValueError(
                "At least one of either 'entries' or 'ref_elts' arguments "
                "have to be given and not empty."
            )

        case (True, False): # 'entries' is given, 'ref_elts' is empty or None
            # Infer reference elements from entries
            # => data for a diagram containing all entries.
            ref_elts = get_elements_from_entries(entries)

        case (False, True): # 'ref_elts' is given, 'entries' is empty or None
            # Init default unary entries for each element
            # => data for a 'blank' diagram with only 0.0 eV entries.
            ref_elts = get_elements(ref_elts)
            entries = get_lacking_elts_entries(entries=[], ref_elts=ref_elts)

        case (True, True): # Both 'entries' and 'ref_elts' are given
            # Filter out entries with unmatching elements
            # Add default unary entries for elements that
            # do not have matching unary entries.
            ref_elts = get_elements(ref_elts)
            entries  = _get_relevant_entries(entries, ref_elts)
            auto_unaries = get_lacking_elts_entries(entries, ref_elts)
            entries += auto_unaries

    return entries, ref_elts


########################################


def _phase_diagram_init(
    entries: Sequence[PDEntry],
    ref_elts: Sequence[Element],
    verbose: bool = False
):
    """Compute a new PhaseDiagram object from given entries and elements."""
    if verbose:
        pd_name = "-".join(list(map(str, ref_elts)))
        print(f"Initializing phase diagram '{pd_name}'")

        new_pd = PhaseDiagram(entries=entries, elements=ref_elts)

        print(f"{pd_name} diagram contains following entries:")
        for entry in new_pd.qhull_entries:
            entry = PDEntry(entry.composition, entry.energy)
            print(f"{entry}, energy_per_atom = {entry.energy_per_atom}")

    else:
        new_pd = PhaseDiagram(entries=entries, elements=ref_elts)

    return new_pd


########################################


def _compute_e_above_hull(
    entries_to_compute: List[PDEntry],
    ref_entries: List[PDEntry],
    stable_limit: float = 0.1,
    verbose: bool = False
) -> List[Dict[str, str|float]]:
    """
    Initialize a phase diagram and compute above hull energies of given entries.
    """
    results = []
    # entries_to_compute is one composition group here,
    # all structures contain the same elements
    comp_refs = get_sub_entries(main_entry=entries_to_compute[0], entry_pool=ref_entries)
    pd_entries, pd_elts = _process_pd_data(entries=comp_refs)
    pd = _phase_diagram_init(pd_entries, pd_elts, verbose=verbose)

    # Compute energy above hull for each generated entry
    for entry in entries_to_compute:
        e_above_hull = pd.get_e_above_hull(entry, allow_negative=True)
        is_stable = e_above_hull <= stable_limit
        if verbose:
            print(f"Entry '{entry.name}':")
            print(f"- energy_per_atom: {entry.energy_per_atom}")
            print(f"- e_above_hull: {e_above_hull}")
            print(f"- is_stable: {is_stable}")
        results.append(
            {
                "name": entry.name,
                "e_above_hull": e_above_hull,
                "stable": is_stable.item()
            }
        )
    return results
#---------------------------------------
def batch_compute_e_above_hull(
    entries: List[List[PDEntry]],
    ref_entries: List[PDEntry],
    stable_limit: float = 0.1,
    workers: int = 1,
    verbose: bool = False
) -> List[Dict[str, Union[str, float]]]:
    """
    For each sublist, build the minimal phase diagram using the reference entries,
    then computes the energy above hull of all entries in the sublist.
    Can be parallelized over phase diagram initialization, one diagram per process.

    Parameters:
        entries ([[PDEntry]]):      List of entry sublists whose energy is needed.
                                    One phase diagram is built for each sublist.

        ref_entries ([PDEntry]):    reference entries used to build the diagrams.
                                    Each diagram is built with only necessary
                                    entries from this sequence.

        stable_limit (float):       Threshold of the energy above hull above which
                                    the entry is considered unstable, in eV/atom.
                                    Defaults to 0.1 eV/atom.

        workers (int):              Number of parallel processes to spawn.

        verbose (bool):             Whether to print each reference entry used when
                                    building a phase diagram.

    Returns:
        List[Dict[str, str|float]]: List of dicts containing the name, energy above hull
        and stringified (for JSON serailization) stability test boolean for one entry
        structure each.
    """
    check_type(entries, "entries", (Sequence,))
    for idx, entry_group in enumerate(entries):
        check_type(entry_group, f"entries[{idx}]", (Sequence,))
    for idx, entry in enumerate(flatten(entries)):
        check_type(entry, f"entries[{idx}]", (PDEntry,))
    check_type(stable_limit, "stable_limit", (float,))
    check_num_value(stable_limit, "stable_limit", ">=", 0.0)
    check_type(workers, "workers", (int,))
    check_num_value(workers, "workers", ">", 0)
    check_type(verbose, "verbose", (bool,))

    energy_computer = partial(
        _compute_e_above_hull,
        ref_entries=ref_entries,
        stable_limit=stable_limit,
        verbose=verbose
    )

    computed_energies = process_map(
        energy_computer,
        entries,
        max_workers=workers,
        chunksize=1,
        desc="Compute above hull energies"
    )
    computed_energies = flatten(computed_energies)
    #computed_energies = []
    #for entry_group in entries:
    #    pd_name = '-'.join([str(elt) for elt in entry_group[0].elements])
    #    print(f"Begin computing phase diagram {pd_name}")
    #    computed_energies += energy_computer(entry_group)
    #    print("Energy computation step complete")

    return computed_energies


########################################
