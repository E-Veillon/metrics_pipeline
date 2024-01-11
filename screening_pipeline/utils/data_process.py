'''
Functions to process raw calculation results from VASP and test actual elimination criterions.
'''


########################################
# SYSTEM I/O MODULES

from typing import Dict, Union, Sequence, Tuple, Any, List, Iterable, Optional
from pathlib import Path

########################################
# OPTIMIZATION MODULES

import numpy as np
from itertools import combinations
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.composition import Composition, Element
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

########################################
# LOCAL MODULES

from screening_pipeline.utils.matcher import group_by_stoichiometry, flatten
from screening_pipeline.utils.periodic_table import get_elements, get_elemental_subsets, \
                                                    get_delta_sol_el_ratio

########################################
# TYPE ALIASES

PathLike    = Union[str, Path]
FormulaLike = Union[str, Iterable[Union[str, int, Element]]]

########################################
# LOCAL FUNCTIONS


# Functions related to phase stability screening with Convex Hull construction method.

def get_elements_from_entries(entries: Sequence[PDEntry]) -> List[Element]:
    '''
    Get a list of all unique elements present in a sequence of PDEntry objects.

    Parameters:
        entries ([PDEntry]): The entries to get unique elements from.

    Returns:
        The list of found Element objects.
    '''
    assert isinstance(entries, Sequence)
    if not entries: return []
    assert all(isinstance(entry, PDEntry) for entry in entries)
    if len(entries) == 1: return entries[0].elements
    return list(set(flatten([entry.elements for entry in entries])))

########################################

def init_entries_and_group_by_dim(
        structs_data: dict[str, dict[str, Any]]
    ) -> List[List[Union[PDEntry, List[PDEntry]]]]:
    '''
    Initialize entries from given data and group them inside nested lists 
    according to the minimal phase diagram necessary for each one.
    Each index of the bigger list represent the dimension 
    (ie. number of distinct elements) of structures inside each sublist.
    Moreover, each sublist contains one subsublist for each composition type.
    In other words, one gets something of the form:
    groups = [
              [], 
              [Elements], 
              [[Binary 1 (eg. all "Fe-O")], [Binary 2 (eg. all "Mn-O")], ...], 
              [[Ternary 1 (eg. all "Fe-Mn-O")], [Ternary 2 (eg. all "Fe-Co-O")], ...], 
              ...
            ]
    
    Parameters:
        structs_data (dict):    structures data as extracted from VASP with 
                                batch_extract_vasp_data function.
    
    Returns:
        All structures data grouped by structure dimensionality and composition.
    '''

    groups = [[]] # fill the index 0 to match indexes and entries dimensionality

    entry_list = [
        PDEntry(
            composition=data['composition'], 
            energy=data['final_energy'], 
            name=name
        ) for name, data in structs_data.items()
    ]

    max_dim  = max([len(entry.elements) for entry in entry_list])
    elements = get_elements_from_entries(entry_list)

    elts_entries = [
        PDEntry(
            composition=Composition(str(elt)), 
            energy=0.
        ) for elt in elements
    ]

    groups.append(elts_entries)
    # At this point, groups = [[], [Elements]]

    for dim in range(2, max_dim + 1):
        group = list(filter(lambda entry: len(entry.elements) == dim, entry_list))
        groups.append(group)
    # At this point,  groups = [[], [Elements], [Binaries], [Ternaries], ...]

    for grp_idx, group in enumerate(groups[2:], start=2):
        groups[grp_idx] = group_by_stoichiometry(group)

    return groups

########################################

def get_sub_entries(
        main_entry: PDEntry, 
        entry_pool: Sequence[PDEntry]
    ) -> List[PDEntry]:
    '''
    Extract all PDEntry objects whose elemental composition are subsets of the
    composition of the main PDEntry object from a sequence of entries.

    Parameters:
        main_entry (PDEntry): The reference entry to search sub-entries from.

        entry_pool ([PDEntry]): The sequence of entries to search sub-entries in.
    
    Returns:
        The list of all found sub-entries relative to the main entry.
    '''

    assert isinstance(main_entry, PDEntry)
    assert isinstance(entry_pool, Sequence)
    if not entry_pool: return []
    assert all(isinstance(entry, PDEntry) for entry in entry_pool)

    sub_entries = list(filter(
        lambda entry: all([elt in main_entry.elements for elt in entry.elements]), 
        entry_pool
    ))

    return sub_entries

########################################

def _get_relevant_entries(
        entries: Sequence[PDEntry], 
        ref_elts: Union[Sequence[Element], set[Element]]
        ) -> List[PDEntry]:
    '''
    Extract entries that are only composed of given reference Elements.

    Parameters:
        entries ([PDEntry]):    Entries to filter out.

        ref_elts ([Elements]):  Reference Elements defining the wanted entry space.
    
    Returns:
        A list of filtered entries.
    '''
    ref_entry = PDEntry(
        composition=Composition([(elt, 1) for elt in ref_elts]), 
        energy=0.0, 
        name='ref_entry'
    )

    return get_sub_entries(main_entry=ref_entry, entry_pool=entries)

########################################

def _auto_define_lacking_elts_entries(
        entries: Sequence[PDEntry], 
        ref_elts: Union[Sequence[Element], set[Element]]
    ) -> List[PDEntry]:
    '''
    Check elemental entries with respect to given reference elements, 
    then auto-initialize lacking elemental entries with an energy of 0.0 eV.

    Parameter:
        entries ([PDEntry]):    Entries to check elemental entries in.

        ref_elts ([Elements]):  Reference Elements defining the wanted entry space.
    
    Returns:
        A list of auto-defined elemental entries.
    '''

    user_elt_entries = list(filter(
    lambda entry: entry.is_element() and entry.elements[0] in ref_elts, 
    entries
    ))

    lacking_elts = list(filter(
        lambda elt: elt not in get_elements_from_entries(user_elt_entries), 
        ref_elts
    ))

    auto_defined_elts_entries = [
        PDEntry(
            composition=Composition(str(elt)), 
            energy=0.0, 
            name=elt.symbol, 
            attribute='auto_elt_ref'
        ) for elt in lacking_elts
    ]

    return auto_defined_elts_entries

########################################

def phase_diagram_init(
        entries: Sequence[PDEntry], 
        ref_elts: Optional[FormulaLike] = None, 
    ) -> PhaseDiagram:
    '''
    Compute a new PhaseDiagram object from given elements and entries.

    Parameters:
        entries ([PDEntry]):        The entries that will be put into the PhaseDiagram.
                                    If some elemental entries are provided, they will be
                                    used with their energy. If some elemental entries are
                                    lacking with respect to ref_elts, they will be auto-
                                    initialized with energy = 0.0 eV. If some entries have
                                    elements not in referenced in provided ref_elts, they
                                    will be ignored.

        ref_elts (str|Iterable):    The  elemental references of the new phase diagram.
                                    If a single string is provided, it can either 
                                    be a raw formula (eg. 'FePO4') or a composition 
                                    string containing element symbols separated by 
                                    '-' (eg. 'Fe-P-O').
                                    If an iterable is given, it can contain valid 
                                    element symbols, atomic numbers and/or Element 
                                    objects.
                                    If not provided, they are computed from given entries.
    Returns:
        The constructed PhaseDiagram object.
    '''

    assert entries and isinstance(entries, Sequence)
    assert all(isinstance(entry, PDEntry) for entry in entries)

    if not ref_elts: ref_elts = get_elements_from_entries(entries)

    else:
        assert isinstance(ref_elts, Iterable)
        ref_elts  = get_elements(ref_elts)
        entries   = _get_relevant_entries(entries, ref_elts)

    entry_list = _auto_define_lacking_elts_entries(entries, ref_elts) + list(entries)

    new_pd = PhaseDiagram(
        entries=entry_list, 
        elements=ref_elts
    )

    return new_pd

########################################
#TODO: finish flexibility update
def _calculate_instability_energies(
        main_entries: Sequence[PDEntry], 
        sub_entries_pool: Sequence[PDEntry]
    ) -> List[Tuple[str, float]]:
    '''
    Initialize a PhaseDiagram from given entries, then calculate relative
    instability energies ΔH for each non-elemental entry in the diagram.

    '''

    energies = []
    ref_elts = get_elements_from_entries(main_entries)

    sub_entries = get_sub_entries(
        main_entry=main_entries[0], 
        entry_pool=sub_entries_pool
    )

    entry_list  = main_entries + sub_entries

    convex_hull = phase_diagram_init(
        entries=entry_list, 
        ref_elts=ref_elts
    )

    for entry in convex_hull.entries[len(ref_elts):]: # skip elemental entries
        delta_H = convex_hull.get_e_above_hull(entry, allow_negative=True)
        energies.append((entry.name, delta_H))

    return energies

########################################
#TODO: finish flexibility update
def calculate_instability_energies(
        entries: Sequence[PDEntry], 
        ref_elts: Optional[Sequence[PDEntry]] = None
    ) -> List[Tuple[str, float]]:
    '''
    Initialize a PhaseDiagram from given entries, then calculate relative
    instability energies ΔH for each non-elemental entry in the diagram.
    The PhaseDiagram's reference elements are automatically initialized 
    from given entries at energy = 0.0 if not provided.

    parameters:
        entries ([PDEntry]):    The entries that will be put into the PhaseDiagram.
                                
        ref_elts ([PDEntry]):   Elemental entries that are references for the PhaseDiagram.
                                If not provided or incomplete, lacking references are
                                automatically initialized with energy = 0.0 eV from given
                                entries. If some provided elements are not used in the entries, 
                                they will be ignored to get minimal diagram dimension and save
                                calculation time and memory.

    Returns:
        List[Tuple[str,float]]: A list of tuples each containing the name of the entry and
                                corresponding ΔH energy in eV/atom.
    '''

    assert isinstance(entries, Sequence)
    if not entries: return []
    assert all(isinstance(entry, PDEntry) for entry in entries)
    
    if ref_elts:
        assert isinstance(ref_elts, Sequence)
        assert all(isinstance(entry, PDEntry) and entry.is_element() for entry in ref_elts)
    
    else: ref_elts = []

    return _calculate_instability_energies(entries, ref_elts)

########################################
#TODO: finish flexibility update
def batch_calculate_instability_energies(structs_data: dict, workers: int = 1):
    '''
    Construct an adaptive convex hull for each structure according to their composition.
    A binary structure does not need comparison with higher order structures.
    However, for a higher order structure, smaller convex hulls can be useful to determine
    its critical formation energy. Therefore, this function constructs the minimal convex hull
    for each compositional group.
    Parameters:
        structs_data (dict):    A dict containing following data about each structure:
                                    - structure directory name (dict's keys), 
                                    - the Structure object, 
                                    - its composition, 
                                    - its relaxed energy (in eV).
    Returns:
        Dict: The same data with all ΔH calculated in 'delta_H' keys.
    '''
    
    assert isinstance(structs_data, dict) and len(structs_data) > 0, \
    f'''Invalid input provided, it either was not a dict or was empty.
        Detected type: {type(structs_data)}.
        Detected length: {len(structs_data)}.'''

    groups              = init_entries_and_group_by_dim(structs_data)
    entry_pool          = []

    for dim_group in groups[2:]:
        if dim_group == []: continue


        setup_calc_inst_energs = partial(
            _calculate_instability_energies, 
            entry_pool=entry_pool
        )

        energies = flatten(list(process_map(
            setup_calc_inst_energs, 
            dim_group, 
            max_workers=workers, 
            chunksize=1, 
            desc=f'computing ΔH for structs of order {groups.index(dim_group)}'
        )))

        for energy in energies:
            structs_data[energy[0]]['delta_H'] = energy[1]

        entry_pool += flatten(dim_group)
    
    return structs_data

########################################

# Functions related to Band Gap screening with Δ-Sol method.

def calculate_delta_sol_band_gap(
        data: dict
    ) -> Tuple[str, float]:
    '''
    Calculate Δ-Sol band gap value of a structure, provided a dict containing all necessary data.

    Parameters:
        data (dict): A dict containing following data about a structure:
                        - Its name, 
                        - The Structure object, 
                        - Its original relaxed total energy, 
                        - Its total energy with more charge density, 
                        - Its total energy with less charge density
    Returns:
        Tuple[str, float]: The name of the structure and its band gap value.
    '''

    assert isinstance(data, dict)

    name         = data['name']
    structure    = data['structure']
    E_N0         = data['E_N0']
    E_N0_plus_n  = data['E_N0_plus_n']
    E_N0_minus_n = data['E_N0_minus_n']
    n_ratio      = get_delta_sol_el_ratio(structure)
    # E_FG = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n -> Δ-Sol band gap 
    # (Ref 32 in screening_pipeline/Bibliography))
    E_band_gap   = (E_N0_plus_n + E_N0_minus_n - 2*E_N0)/n_ratio

    return name, E_band_gap

########################################

def batch_calculate_delta_sol_band_gaps(
        structs_data: dict, 
        bg_structs_data: dict, 
        workers: int = 1, 
        /
    ) -> Dict[str, float]:
    '''
    Calculate Δ-Sol band gap value for every structure in a batch from their data,
    as provided by extract_vasp_data_for_delta_sol function applied on the 3 energy calculations.

    Parameters:
        structs_data (dict):    Dict containing the original data extracted from previous VASP calculation.

        bg_structs_data (dict): Dict of the calculated data on structures with changed charge density.

        workers (int):          The number of parallel processes to spawn. Defaults to 1.
    
    Returns:
        Dict[str, float]: Dict of Band gap values associated with the original structure directory name.
    '''

    assert isinstance(structs_data, dict)
    assert isinstance(bg_structs_data, dict)
    assert isinstance(workers, int) and workers >= 0

    final_energies = {
        name: {
            'name': name, 
            'structure': data['structure'], 
            'E_N0': data['final_energy'], 
            'E_N0_plus_n': bg_structs_data['_'.join(name, 'plus')]['final_energy'], 
            'E_N0_minus_n': bg_structs_data['_'.join(name, 'minus')]['final_energy']
        } for name, data in structs_data.items()
    }

    nbr_structs = len(final_energies)
    chunksize   = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)
    data_list   = list(final_energies.values())

    E_band_gaps = list(process_map(
        calculate_delta_sol_band_gap, 
        data_list, 
        workers=workers, 
        chhunksize=chunksize
    ))

    E_band_gaps = dict(E_band_gaps)

    return E_band_gaps

########################################

def filter_by_band_gap(
        E_band_gaps: dict, 
        valid_interval: Sequence[float], 
        base_dir: PathLike, 
        ignore_file: str
    ) -> None:

    for name, E_band_gap in E_band_gaps.items():

        name_plus    = '_'.join(name, 'plus')
        name_minus   = '_'.join(name, 'minus')
        bg_too_small = E_band_gap < min(valid_interval)
        bg_too_big   = E_band_gap > max(valid_interval)

        if bg_too_small or bg_too_big:

            reject_str = f'Δ-Sol band gap was estimated to {E_band_gap} eV, \
                        which is not inside the interval [{min(valid_interval)}, {max(valid_interval)}].\n \
                        Therefore, it is not suitable for wanted application, \
                        it should not be considered in further screening steps.'
        
            path_plus  = Path(base_dir / name_plus / ignore_file)
            path_plus.touch()
            path_plus.write_text(reject_str)
        
            path_minus = Path(base_dir / name_minus / ignore_file)
            path_minus.touch()
            path_minus.write_text(reject_str)

########################################