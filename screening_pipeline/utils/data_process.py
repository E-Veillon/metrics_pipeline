"""
Functions to process raw calculation results from VASP and test actual elimination criterions.
"""

########################################
# SYSTEM I/O MODULES

from typing import Dict, Union, Sequence, Tuple, Any, List, Iterable, Optional, Literal
from pathlib import Path

########################################
# OPTIMIZATION MODULES

import numpy as np
from itertools import combinations
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.core.structure import SiteCollection
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

########################################
# LOCAL MODULES

from screening_pipeline.utils.utils import flatten
from screening_pipeline.utils.custom_types import PathLike, FormulaLike
from screening_pipeline.utils.matcher import group_by_composition
from screening_pipeline.utils.periodic_table import get_elements, get_delta_sol_el_ratio

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


def get_max_dim(entries: Sequence[PDEntry]) -> int:
    """Search for the maximum number of distinct elements in given entries."""
    return max([len(entry.elements) for entry in entries])


########################################


def init_entries_from_dict(
        entries_dict: Dict[str, Dict], attribute: Optional[str] = None, workers: int = 1
    ) -> List[PDEntry]:
    """
    Convert structures data dict into a list of PDEntry objects compatible with phase diagrams.

    Parameters:
        entries_dict (dict):    Dict of dicts containing following structure data keys:
                                - entry_id (either dataset ID or generated name),
                                - composition (Composition object),
                                - final_energy (float).

        attribute (str):        Label to put as attribute in the entries.
    
    Returns:
        List[PDEntry]: The list of phase diagram entries.
    """
    if not isinstance(entries_dict, Dict):
        raise TypeError(
            f"'entries_dict' arg expected a type 'dict', got '{type(entries_dict)}' instead."
        )
    
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
        entries: Dict[str,Any]|List[PDEntry],
        max_dim: Optional[int] = None,
        ref_elts: Optional[List[Element]] = None,
        workers: int = 1
    ) -> List[PDEntry]:
    """
    Filter out entries in database file that are useless for generated data.
    Data format conversion can be parallelized.

    Parameters:
        entries (dict):         Entries to filter.

        max_dim (int):          Max dimension of entries to keep.

        ref_elts ([Element]):   List of useful Elements to keep in the dataset.
                                All structures containing other elements are discarded.
        
        workers (int):          Number of parallel processes to spawn.

    Returns:
        List[PDEntry]: List of useful entries.
    """

    assert isinstance(entries, (Dict, List))
    assert isinstance(max_dim, int) or max_dim is None
    assert isinstance(ref_elts, List) or ref_elts is None

    print("Filtering reference dataset...")

    if isinstance(entries, Dict):
        entries = init_entries_from_dict(entries, attribute="ref_struct", workers=workers)

    entries = _get_relevant_entries(entries, ref_elts)
    entries = list(filter(
        lambda entry: len(entry.elements) <= max_dim,
        entries
    ))
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

    assert isinstance(entries, Sequence)
    if not entries:
        return []
    assert all(isinstance(entry, PDEntry) for entry in entries)
    if len(entries) == 1:
        return entries[0].elements

    return list(set(flatten([entry.elements for entry in entries])))


########################################


def group_by_dim_and_comp(
        entries: Sequence[PDEntry], max_dim: Optional[int] = None, workers: int = 1
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

    def _group_one_dim(entries: Sequence[PDEntry], dim: int):
        dim_group = list(filter(lambda entry: len(entry.elements) == dim, entries))
        grouped_entries = group_by_composition(dim_group)
        return grouped_entries

    initial_group = [[]]  # fill the index 0 to match indexes and entries dimensionality

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


#def init_entries_and_group_by_dim_and_comp(
#    structs_data: Dict[str, Dict[str, Any]],
#    ref_structs: Optional[Dict[str, Dict[str, Any]]] = None,
#) -> List[List[Union[PDEntry, List[PDEntry]]]]:
#    """
#    Initialize entries from given data and group them inside nested lists
#    according to the minimal phase diagram necessary for each one.
#    Each index of the bigger list represent the dimension
#    (ie. number of distinct elements) of structures inside each sublist.
#    Moreover, each sublist contains one subsublist for each composition type.
#    In other words, one gets something of the form:
#    groups = [
#              [],
#              [[Unary 1 (e.g. all "Fe")], [Unary 2 (e.g. all "Mn")]], ...],
#              [[Binary 1 (e.g. all "Fe-O")], [Binary 2 (e.g. all "Mn-O")], ...],
#              [[Ternary 1 (e.g. all "Fe-Mn-O")], [Ternary 2 (e.g. all "Fe-Co-O")], ...],
#              ...etc
#            ]
#
#    Parameters:
#        structs_data (dict):    Structure data as extracted from VASP with
#                                batch_extract_vasp_data function. This arg
#                                is for relaxed structures that need a ΔH
#                                computation (generated structures).
#
#        ref_structs (dict):     Structure data as extracted from VASP with
#                                batch_extract_vasp_data function. This arg
#                                is for structures extracted from a dataset,
#                                that can give some reference about previously
#                                found energies in the space of interest.
#                                The ΔH energy of these structures will not
#                                be computed nor returned. Defaults to None.
#
#    Returns:
#        All structures data grouped by structure dimensionality and composition.
#    """
#
#    if not isinstance(structs_data, dict):
#        raise TypeError(
#            "'struct_data' argument expected a 'dict', "
#            f"got '{type(structs_data)}' instead."
#        )
#    if not structs_data:
#        return []
#    assert all(
#        isinstance(name, str) and isinstance(data, dict)
#        for name, data in structs_data.items()
#    )
#    print("Converting generated data to entries...")
#    entry_list = [
#        PDEntry(
#            composition=data["composition"],
#            energy=data["first_ionic_energy"],
#            name=name,
#            attribute="generated",
#        )
#        for name, data in structs_data.items()
#    ]
#    print("Conversion done.")
#
#    max_dim = get_max_dim(entry_list)
#    elements = get_elements_from_entries(entry_list)
#
#    if ref_structs:
#        assert isinstance(ref_structs, dict)
#        assert all(
#            isinstance(name, str) and isinstance(data, dict)
#            for name, data in ref_structs.items()
#        )
#        print("reference dataset detected.")
#
#        ref_entry_list = filter_database_entries(
#            entries=ref_structs, max_dim=max_dim, ref_elts=elements
#        )
#    else:
#        ref_entry_list = []
#
#    elt_entries = get_lacking_elts_entries(
#        entries=entry_list, ref_elts=elements
#    )
#
#    full_entry_list = entry_list + ref_entry_list + elt_entries
#
#    print("Grouping data by dim and composition...")
#    groups = [[]]  # fill the index 0 to match indexes and entries dimensionality
#
#    for dim in range(1, max_dim + 1):
#        print(f"filling dim group {dim}...")
#        group = list(filter(lambda entry: len(entry.elements) == dim, full_entry_list))
#        group = group_by_composition(group)
#        groups.append(group)
#        print("group filled.")
#    print("Grouping done.")
#
#    return groups


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

    assert isinstance(main_entry, PDEntry)
    assert isinstance(entry_pool, Sequence)
    if not entry_pool:
        return []
    assert all(isinstance(entry, PDEntry) for entry in entry_pool), (
        "Some of given entries in 'entry_pool' arguments are not instances of PDEntry class.\n"
        "See below the set of types found in the argument:\n"
        f"{set([type(entry) for entry in entry_pool])}"
    )

    sub_entries = list(
        filter(
            lambda entry: all([elt in main_entry.elements for elt in entry.elements]),
            entry_pool,
        )
    )

    return sub_entries


########################################


def _get_relevant_entries(
    entries: Sequence[PDEntry], ref_elts: Union[Sequence[Element], set[Element]]
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
    entries: Sequence[PDEntry], ref_elts: Union[Sequence[Element], set[Element]]
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

    if not entries:
        user_elt_entries = []
    else:
        user_elt_entries = list(
            filter(
                lambda entry: entry.is_element and entry.elements[0] in ref_elts,
                entries
            )
        )

    lacking_elts = list(
        filter(
            lambda elt: elt not in get_elements_from_entries(user_elt_entries),
            ref_elts
        )
    )

    auto_defined_elts_entries = [
        PDEntry(
            composition=Composition(str(elt)),
            energy=0.0,
            name=elt.symbol,
            attribute="element_ref",
        )
        for elt in lacking_elts
    ]

    return auto_defined_elts_entries


########################################


def phase_diagram_init(
    entries: Sequence[PDEntry],
    ref_elts: Optional[FormulaLike] = None,
) -> PhaseDiagram:
    """
    Compute a new PhaseDiagram object from given elements and entries.

    Parameters:
        entries ([PDEntry]):        The entries that will be put into the PhaseDiagram.
                                    If some elemental entries are provided, they will be
                                    used with their energy. If some elemental entries are
                                    lacking with respect to ref_elts, they will be
                                    initialized with energy = 0.0 eV. If some entries have
                                    elements not referenced in provided ref_elts, they
                                    will be ignored. This behaviour is particularly
                                    useful if one needs to initialize several diagrams from
                                    different parts of the same entry dataset.

        ref_elts (str|Sequence):    The  elemental references of the new phase diagram.
                                    If a single string is provided, it can either
                                    be a raw formula (eg. 'FePO4') or a composition
                                    string containing element symbols separated by
                                    '-' (eg. 'Fe-P-O').
                                    If a sequence is given, it can contain valid
                                    element symbols, atomic numbers and/or Element
                                    objects.
                                    If not provided, they are computed from given entries.
                                    In that case, all provided entries are checked, so the
                                    phase diagram will be of the minimal dimension that
                                    contains all entries.
    Returns:
        The constructed PhaseDiagram object.
    """

    assert entries and isinstance(entries, Sequence)
    assert all(isinstance(entry, PDEntry) for entry in entries)

    if not ref_elts:
        ref_elts = get_elements_from_entries(entries)

    else:
        assert isinstance(ref_elts, Sequence)
        ref_elts = get_elements(ref_elts)
        entries = _get_relevant_entries(entries, ref_elts)

    entry_list = get_lacking_elts_entries(entries, ref_elts) + list(entries)

    pd_name = "-".join(list(map(str, ref_elts)))
    print(f"Initializing phase diagram '{pd_name}'")
    new_pd = PhaseDiagram(entries=entry_list, elements=ref_elts)
    print(f"{pd_name} diagram contains following entries:")
    for entry in new_pd.qhull_entries:
        print(f"{entry}")

    return new_pd


########################################


def _compute_e_above_hull(
        entries_to_compute: List[PDEntry], ref_entries: List[PDEntry],
        stable_limit: float = 0.1
    ) -> List[Dict[str, str|float]]:
    """
    Initialize a phase diagram and compute above hull energies of given entries.
    """
    results = []
    comp_refs = get_sub_entries(main_entry=entries_to_compute[0], entry_pool=ref_entries)
    pd = phase_diagram_init(entries=comp_refs)

    # Compute energy above hull for each generated entry
    for entry in entries_to_compute:
        e_above_hull = pd.get_e_above_hull(entry, allow_negative=True)
        is_stable = e_above_hull <= stable_limit
        print(f"Entry '{entry.name}':")
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
        entries_to_compute: List[List[PDEntry]],
        ref_entries: List[PDEntry],
        stable_limit: float = 0.1,
        workers: int = 1
    ) -> List[Dict[str, Union[str, float]]]:
    """
    For each sublist, build the minimal phase diagram using the reference entries,
    then computes the energy above hull of all entries in the sublist.
    Can be parallelized over phase diagram initialization, one diagram per process.

    Parameters:
        entries_to_compute ([[PDEntry]]):   List of entry sublists whose energy is needed.
                                            One phase diagram is built for each sublist.

        ref_entries ([PDEntry]):            reference entries used to build the diagrams.
                                            Each diagram is built with only necessary
                                            entries from this sequence.

        stable_limit (float):               Threshold of the energy above hull above which
                                            the entry is considered unstable, in eV/atom.
                                            Defaults to 0.1 eV/atom.

        workers (int):                      Number of parallel processes to spawn.
        
    Returns:
        List[Dict[str, str|float]]: List of dicts containing the name, energy above hull
        and stringified (for JSON serailization) stability test boolean for one entry
        structure each.
    """

    energy_computer = partial(
        _compute_e_above_hull,
        ref_entries=ref_entries,
        stable_limit=stable_limit
    )

    computed_energies = process_map(
        energy_computer,
        entries_to_compute,
        max_workers=workers,
        chunksize=1,
        desc="Compute above hull energies"
    )
    computed_energies = flatten(computed_energies)

    return computed_energies


########################################


#def _calculate_instability_energies(
#    main_entries: Sequence[PDEntry], sub_entries_pool: Sequence[PDEntry]
#) -> List[Tuple[str, float]]:
#    """
#    Initialize a PhaseDiagram from given entries, then calculate relative
#    instability energies ΔH for each generated entry in the diagram, ie.
#    non-elemental nor reference structures. Entries that need to have
#    their ΔH calculated must have the string 'generated' as entry.attribute.
#    The PhaseDiagram's reference elements are automatically initialized
#    from given entries at energy = 0.0 eV if not provided.
#
#    Parameters:
#        main_entries ([PDEntry]):       Entries of the highest dimension that will
#                                        serve as reference to define the space of
#                                        the diagram.
#
#        sub_entries_pool ([PDEntry]):   All other entries that should be put into
#                                        the diagram. If an entry in this argument
#                                        contains elements that are not referenced
#                                        in any main entry, it will not be put in
#                                        the diagram and its energy will not be
#                                        computed.
#
#    Returns:
#        List[Tuple[str,float]]: A list of tuples each containing the name of the entry
#                                and corresponding ΔH energy in eV/atom.
#    """
#
#    energies = []
#    ref_elts = get_elements_from_entries(main_entries)
#
#    sub_entries = _get_relevant_entries(entries=sub_entries_pool, ref_elts=ref_elts)
#
#    entry_list = list(set(main_entries + sub_entries))
#
#    diagram_entries = list(
#        filter(lambda entry: entry.attribute != "generated", entry_list)
#    )
#    generated_entries = list(
#        filter(lambda entry: entry.attribute == "generated", entry_list)
#    )
#    if not generated_entries:
#        comp = "-".join(list(map(str, ref_elts)))
#        return energies
#
#    convex_hull = phase_diagram_init(entries=diagram_entries, ref_elts=ref_elts)
#
#    for entry in generated_entries:
#        delta_H = convex_hull.get_e_above_hull(entry, allow_negative=True)
#        energies.append((entry.name, delta_H))
#
#    return energies


########################################


#def calculate_instability_energies(
#    entries: Sequence[PDEntry], ref_elts: Optional[Sequence[PDEntry]] = None
#) -> List[Tuple[str, float]]:
#    """
#    Initialize a PhaseDiagram from given entries, then calculate relative
#    instability energies ΔH for each generated entry in the diagram, ie.
#    non-elemental nor reference structures. Entries that need to have
#    their ΔH calculated must have the string 'generated' as entry.attribute.
#    The PhaseDiagram's reference elements are automatically initialized
#    from given entries at energy = 0.0 eV if not provided.
#
#    parameters:
#        entries ([PDEntry]):    The entries that will be put into the PhaseDiagram.
#
#        ref_elts ([PDEntry]):   Elemental entries that are references for the PhaseDiagram.
#                                If not provided or incomplete, lacking references are
#                                automatically initialized with energy = 0.0 eV from given
#                                entries. If some provided elements are not used in the entries,
#                                they will be ignored to get minimal diagram dimension and save
#                                calculation time and memory.
#                                This argument is specifically designed in case non-zero elemental
#                                energy references are needed, and can be ignored in other cases.
#
#    Returns:
#        List[Tuple[str,float]]: A list of tuples each containing the name of the entry and
#                                corresponding ΔH energy in eV/atom.
#    """
#
#    assert isinstance(entries, Sequence)
#    if not entries:
#        return []
#    assert all(isinstance(entry, PDEntry) for entry in entries)
#
#    if ref_elts:
#        assert isinstance(ref_elts, Sequence)
#        assert all(
#            isinstance(entry, PDEntry) and entry.is_element for entry in ref_elts
#        )
#
#    else:
#        ref_elts = []
#
#    return _calculate_instability_energies(entries, ref_elts)


########################################


#def batch_calculate_instability_energies(
#    structs_data: dict,
#    structs_ref: Optional[Dict[str, Dict[str, Any]]] = None,
#    workers: int = 1,
#):
#    """
#    Construct an adaptive convex hull for each structure according to their composition.
#    A binary structure does not need comparison with higher order structures.
#    However, for a higher order structure, convex hulls of smaller order can be useful to
#    determine its critical formation energy. Therefore, this function constructs the minimal
#    convex hull for each compositional group.
#
#    Parameters:
#        structs_data (dict):    A dict containing following data about each structure:
#                                    - its name (dict's keys),
#                                    - its composition as a Composition object,
#                                    - its energy after the first ionic step,
#                                    - its relaxed energy (in eV).
#
#        structs_ref (dict):     Reference data as a dictionnary where each key is the name
#                                of the structure and each value is another dictionnary containing
#                                the composition and total energy.
#
#    Returns:
#        Dict: The same structs_data dict with all ΔH calculated in 'delta_H' keys.
#    """
#
#    assert isinstance(structs_data, dict) and len(structs_data) > 0, (
#        "Invalid input provided, it either was not a dict or was empty.\n"
#        f"Detected type: {type(structs_data)}.\n"
#        f"Detected length: {len(structs_data)}.\n"
#    )
#    dim_groups = init_entries_and_group_by_dim_and_comp(
#        structs_data, ref_structs=structs_ref
#    )
#
#    entry_pool = dim_groups[0]
#
#    for dim_group in dim_groups[1:]:
#        if dim_group == []:
#            continue
#
#        setup_calc_inst_energs = partial(
#            _calculate_instability_energies, sub_entries_pool=entry_pool
#        )
#
#        energies = flatten(
#            list(
#                process_map(
#                    setup_calc_inst_energs,
#                    dim_group,
#                    max_workers=workers,
#                    chunksize=1,
#                    desc=f"computing ΔH for structs of order {dim_groups.index(dim_group)}",
#                )
#            )
#        )
#
#        for energy in energies:
#            name = energy[0]
#            delta_H = energy[1]
#            structs_data[name]["delta_H"] = delta_H
#
#        entry_pool += flatten(dim_group)
#
#    return structs_data


########################################
# Functions related to Band Gap screening with Δ-Sol method.


def calculate_delta_sol_band_gap(data: dict) -> Union[Tuple[str, float], Tuple[str, float, float, float]]:
    """
    Calculate Δ-Sol band gap value of a structure, provided a dict containing all necessary data.

    Parameters:
        data (dict): A dict containing at least following data about a structure:
                        - Its name,
                        - The DFT functional used for the calculation, 
                        - The Structure object, 
                        - Its total energy with N0 electrons, 
                        - Its total energy with N0 + n(best) electrons, 
                        - Its total energy with N0 - n(best) electrons.

                     It can also contain data for uncertainty calculations:
                        - Total energy with N0 + n(min) electrons,
                        - Total energy with N0 - n(min) electrons,
                        - Total energy with N0 + n(max) electrons,
                        - Total energy with N0 - n(max) electrons.

    Returns:
        The name of the structure and its band gap value(s).
    """

    def has_str_key(dct: Dict, key: str) -> bool:
        return dct.get(key) is not None

    data_keys = ("name","functional","structure","E_N0","E_N0_plus_n_best","E_N0_minus_n_best")
    supp_keys = ("E_N0_plus_n_min","E_N0_minus_n_min","E_N0_plus_n_max","E_N0_minus_n_max")

    assert isinstance(data, dict)
    assert all([has_str_key(data, key) for key in data_keys])

    n_ratio_best = get_delta_sol_el_ratio(data["structure"], dft_functional=data["functional"], n_star_type="BEST")

    # E_FG = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n -> Δ-Sol band gap 
    # (Ref 32 in screening_pipeline/Bibliography))
    E_band_gap = (data["E_N0_plus_n_best"] + data["E_N0_minus_n_best"] - 2*data["E_N0"])/n_ratio_best

    if not all([has_str_key(data, key) for key in supp_keys]): # data do not have uncertainty keys
        return (data["name"], E_band_gap)

    n_ratio_min = get_delta_sol_el_ratio(data["structure"], dft_functional=data["functional"], n_star_type="MIN")
    n_ratio_max = get_delta_sol_el_ratio(data["structure"], dft_functional=data["functional"], n_star_type="MAX")

    E_band_gap_min = (data["E_N0_plus_n_min"] + data["E_N0_minus_n_min"] - 2*data["E_N0"])/n_ratio_min
    E_band_gap_max = (data["E_N0_plus_n_max"] + data["E_N0_minus_n_max"] - 2*data["E_N0"])/n_ratio_max

    return (data["name"], E_band_gap, E_band_gap_min, E_band_gap_max)


########################################


def batch_calculate_delta_sol_band_gaps(
        bg_data: dict, 
        dft_functional: Literal["LDA","PBE","AM05"] = "PBE", 
        with_uncertainties: bool = False, 
        workers: int = 1, 
        /
    ) -> Dict[str, float]:
    """
    Calculate Δ-Sol band gap value for every structure in a batch from their data,
    as provided by extract_vasp_data_for_delta_sol function applied on the 3 energy calculations.

    Parameters:
        bg_data (dict):               Dict containing structures data extracted from previous VASP static calculations.

        dft_functional (str):         The type of functional used for static calculations. Supported functionals are
                                      "LDA", "PBE", and "AM05". Defaults to "PBE".

        with_uncertainties (bool):    Whether to include uncertainty calculations data in the results.
                                      Defaults to False.

        workers (int):                The number of parallel processes to spawn. Defaults to 1.

    Returns:
        Dict[str, float]: Dict of Band gap values associated with the original structure directory name.
    """

    assert isinstance(bg_data, dict)
    assert dft_functional in {"LDA", "PBE", "AM05"}
    assert isinstance(with_uncertainties, bool)
    assert isinstance(workers, int) and workers >= 1

    if not with_uncertainties:
        final_energies = {
            name: {
                "name": name,
                "functional": dft_functional,
                "structure": data["structure"],
                "E_N0": data[name + "_neutral"],
                "E_N0_plus_n_best": data[name + "_best_plus"],
                "E_N0_minus_n_best": data[name + "_best_minus"],
            } for name, data in bg_data.items()
        }
    
    else:
        final_energies = {
            name: {
                "name": name,
                "functional": dft_functional,
                "structure": data["structure"],
                "E_N0": data[name + "_neutral"],
                "E_N0_plus_n_best": data[name + "_best_plus"],
                "E_N0_minus_n_best": data[name + "_best_minus"],
                "E_N0_plus_n_min": data[name + "_min_plus"],
                "E_N0_minus_n_min": data[name + "_min_minus"],
                "E_N0_plus_n_max": data[name + "_max_plus"],
                "E_N0_minus_n_max": data[name + "_max_minus"],
            } for name, data in bg_data.items()
        }

    nbr_structs = len(final_energies)
    chunksize = min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1
    data_list = list(final_energies.values())

    E_band_gaps = list(process_map(
        calculate_delta_sol_band_gap, 
        data_list, 
        max_workers=workers, 
        chunksize=chunksize
    ))

    if not with_uncertainties:
        E_band_gaps = {
            tup[0]: {
                'E_band_gap': tup[1]
            } for tup in E_band_gaps
        }

    else:
        E_band_gaps = {
            tup[0]: {
                'E_band_gap': tup[1], 
                'E_band_gap_min': tup[2], 
                'E_band_gap_max': tup[3]
            } for tup in E_band_gaps
        }

    return E_band_gaps


########################################
