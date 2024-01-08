'''
Functions to process raw calculation results from VASP and test actual elimination criterions.
'''


########################################
# SYSTEM I/O MODULES

from typing import Dict, Union, Sequence, Tuple, Any, List
from pathlib import Path

########################################
# OPTIMIZATION MODULES

from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE


########################################
# LOCAL MODULES

from screening_pipeline.utils.matcher import group_by_stoichiometry, flatten
from screening_pipeline.utils.periodic_table import get_delta_sol_el_ratio

########################################
# TYPE ALIASES

PathLike = Union[str, Path]

########################################
# LOCAL FUNCTIONS

########################################
# stability screening (Convex Hull construction)

def group_by_dim(structs_data: dict[str, dict[str, Any]]) -> List[List[List[Tuple[str, dict[str, Any]]]]]:
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
                                extract_vasp_data_for_convex_hull function.
    
    Returns:
        All structures data grouped by structure dimensionality and composition.
    '''
    groups       = [[]] # fill the index 0 to have correspondance between index and dim
    structs_list = [(name, data) for name, data in structs_data.items()]
    max_elts_nbr = max([len(data['composition']) for data in structs_data.values()])
    elements     = list(set(flatten([struct[1]['composition'].elements for struct in structs_list])))
    groups.append(elements)

    for elts_nbr in range(2, max_elts_nbr + 1):
        group = list(filter(lambda struct: len(struct[1]['composition']) == elts_nbr, structs_list))
        groups.append(group)
    # At this point,  groups = [[], [Elements], [Binaries], [Ternaries], ...]

    for grp_idx, group in enumerate(groups[2:], start=2):
        groups[grp_idx] = group_by_stoichiometry(group)

    return groups

########################################
# Band Gap screening (Δ-Sol method)

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