'''
Functions to process raw calculation results from VASP and test actual elimination criterions.
'''


########################################
# SYSTEM I/O MODULES

from typing import Dict, Union, Sequence, Tuple, Any, List, Iterable
from pathlib import Path

########################################
# OPTIMIZATION MODULES

import numpy as np
from itertools import combinations
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

def group_by_dim(structs_data: dict[str, dict[str, Any]]) -> List[List[Union[PDEntry, List[PDEntry]]]]:
    '''
    Groups structures inside nested lists according to the minimal 
    phase diagram necessary for each one.
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
    elements = list(set(flatten([entry.elements for entry in entry_list])))

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

def phase_diagram_init(
        ref_elts: FormulaLike, 
        entries: Sequence[PDEntry]
    ) -> PhaseDiagram:
    '''
    Compute a new PhaseDiagram object from given elements and structure data.

    Parameters:
        ref_elts (str|Iterable):        The  elemental references of the new phase diagram.
                                        If a single string is provided, it can either 
                                        be a raw formula (eg. 'FePO4') or a composition 
                                        string containing element symbols separated by 
                                        '-' (eg. 'Fe-P-O').
                                        If an iterable is given, it can contain valid 
                                        element symbols, atomic numbers and/or Element 
                                        objects.

        structs_data (dict|Sequence):   The structures data to put into the diagram.
    
    Returns:
        The constructed PhaseDiagram object.
    '''

    assert all(isinstance(entry, PDEntry) for entry in entries)

    ref_elts     = get_elements(ref_elts)
    elts_entries = [
        PDEntry(
            composition=Composition(str(elt)), 
            energy=0.0, 
            name=elt.symbol, 
            attribute='element_ref'
        ) for elt in ref_elts
    ]

    entry_list = elts_entries + list(entries)

    new_pd = PhaseDiagram(
        entries=entry_list, 
        elements=ref_elts
    )

    return new_pd

########################################

def get_sub_entries(
        main_entry: PDEntry, 
        entry_pool: Sequence[PDEntry]
    ) -> List[PDEntry]:
    '''
    Extract all PDEntry objects which elemental composition are subsets of the
    composition of the main PDEntry object from a Sequence of entries.

    Parameters:
        main_entry (PDEntry): The reference entry to search sub-entries from.

        entry_pool ([PDEntry]): The sequence of entries to search sub-entries in.
    
    Returns:
        The list of all found sub-entries relative to the main entry.
    '''

    assert isinstance(main_entry, PDEntry)
    assert all(isinstance(entry, PDEntry) for entry in entry_pool)

    sub_entries = list(filter(
        lambda entry: all([elt in main_entry.elements for elt in entry.elements]), 
        entry_pool
    ))

    return sub_entries

########################################
# TODO: Work in Progress, it does not work yet.
def init_pd_from_cache(
        ref_elts: FormulaLike, 
        cached_pd_data: dict, 
        new_data: dict[str, dict]|Sequence[Tuple]
    ) -> PhaseDiagram:
    '''
    Initialize a higher order PhaseDiagram object by combining data from lesser ones already computed.
    Avoid expensive construction of complex phase diagrams and redundance in computations.
    If some parts of the higher diagram are not computed yet, this function computes them from scratch
    after combining cached data as much as possible to complete the diagram.

    Parameters:
        ref_elts (str|Iterable):        The  elemental references of the new phase diagram.
                                        If a single string is provided, it can either 
                                        be a raw formula (eg. 'FePO4') or a composition 
                                        string containing element symbols separated by 
                                        '-' (eg. 'Fe-P-O').
                                        If an iterable is given, it can contain valid 
                                        element symbols, atomic numbers and/or Element 
                                        objects.
        
        cached_pd_data (dict):          A dict referencing computed diagrams by their name, 
                                        and containing their MSONable dicts, as returned by
                                        PhaseDiagram.as_dict() method.
        
        new_data (dict|[Tuple]):        Structure data to put in the phase diagram in addition to 
                                        previous diagrams, typically structures of same elemental
                                        composition as the new diagram itself. 
                                        They must be provided as a dict of {name: data} or a sequence 
                                        of (name, data) tuples, where 'name' is the name of the 
                                        structure directory where data were took from, and 'data' is 
                                        a dict containing same infos as provided by the function
                                        'screening_pipeline.utils.vasp_io.extract_vasp_data_for_convex_hull'.
    
    Returns:
        The constructed PhaseDiagram object.
    '''

    assert isinstance(ref_elts, (str, Iterable))
    assert isinstance(cached_pd_data, dict)
    assert isinstance(new_data, (dict, Sequence))

    ref_elts        = get_elements(ref_elts)
    new_pd_dim      = len(ref_elts)
#    new_pd_edges    = set(combinations(ref_elts, 2))
#    satisfied_edges = set()


    all_sub_pd_list = get_elemental_subsets(ref_elts, cached_pd_data.keys())
    all_entries = list(set([cached_pd_data[sub_pd].all_entries for sub_pd in all_sub_pd_list]))
    qhull_entries = list(set([cached_pd_data[sub_pd].qhull_entries for sub_pd in all_sub_pd_list]))
    qhull_data = np.array(list(set([cached_pd_data[sub_pd].qhull_data for sub_pd in all_sub_pd_list])))
    facets = list(set([cached_pd_data[sub_pd].facets for sub_pd in all_sub_pd_list]))
    simplexes = list(set([cached_pd_data[sub_pd].simplexes for sub_pd in all_sub_pd_list]))

    new_pd_data  = {
        "@module": PhaseDiagram.__module__, # OK
        "@class": PhaseDiagram.__name__, # OK
        "all_entries": [entry.as_dict() for entry in all_entries], 
        "elements": [elt.as_dict() for elt in ref_elts], # OK
        "computed_data": {
            "facets": facets, # utiliser get_facets une fois tout les points dans le tableau
            "simplexes": simplexes, # transformer les facets en Simplex et le lister ici
            "all_entries": all_entries, # à calculer
            "qhull_data": qhull_data, # numpy.ndarray, à caster en liste et remettre en array ensuite
            "dim": new_pd_dim, # OK
            "el_refs": [ # OK
                (elt, PDEntry(
                    composition=Composition(elt), 
                    energy=0.0, 
                    name=elt.symbol, 
                    attribute='element_ref'
                )) for elt in ref_elts
            ], 
            "qhull_entries": qhull_entries # à calculer
        }
    }

    '''for sub_dim in reversed(range(2, new_pd_dim)):

        dim_sub_pd_list = list(filter(
            lambda sub_pd: len(get_elements(sub_pd)) == sub_dim, 
            all_sub_pd_list
        ))

        if not dim_sub_pd_list:
            continue

        # TODO: Il faut éviter d'ajouter des sous-diagrammes si les arêtes 
        #       correspondantes sont déjà satisfaites.
        for sub_pd in dim_sub_pd_list:
            sub_pd_edges = set(combinations(get_elements(sub_pd), 2))
            for edge in sub_pd_edges:
                if edge in satisfied_edges:
                    dim_sub_pd_list.remove(sub_pd)

        if not dim_sub_pd_list:
            continue

        new_pd_data['all_entries'] += list(set(sum(
            cached_pd_data[pd_name]['all_entries'] for pd_name in dim_sub_pd_list
        )))

    new_pd_data['computed_data']['all_entries'] = [
        PDEntry.from_dict(entry) for entry in new_pd_data['all_entries']
    ]'''

    new_pd = PhaseDiagram.from_dict(dct=new_pd_data)
    
    return new_pd

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