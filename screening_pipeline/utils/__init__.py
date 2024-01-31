from .periodic_table import has_rare_gas, discard_rare_gas_structures, \
    has_rare_earth, discard_rare_earth_structures, get_elements, get_elemental_subsets, \
    get_all_elements_groups, get_element_valence_electrons, get_all_valence_electrons
from .cif_io import read_cif, write_cif, extract_cif_from_file, cif_str_to_struct
from .matcher import remove_equivalent
from .vasp_io import vasp_relaxation_settings, vasp_static_settings, \
                    vasp_launcher, vasp_batch_launch, \
                    batch_extract_vasp_data, chgcar_density_switch, delta_sol_inputs_init
from .spacegroup import structure_symmetrizer, batch_symmetrizer
from .fitted_values import E_O2_FIT, U_VALUES, DELTA_E_M, EXP_DELTA_H, EL_PER_XC_VOL
from .paths import add_new_dir, batch_add_new_dirs
from .data_process import get_elements_from_entries, init_entries_and_group_by_dim_and_comp, \
    get_sub_entries, phase_diagram_init, calculate_instability_energies, \
    batch_calculate_instability_energies, calculate_delta_sol_band_gap, \
    batch_calculate_delta_sol_band_gaps

__all__ = [
    "has_rare_gas", "discard_rare_gas_structures", 
    "has_rare_earth", "discard_rare_earth_structures", 
    "get_elements", "get_elemental_subsets", "get_all_elements_groups", "get_all_valence_electrons", 
    "read_cif", "write_cif", 
    "structure_symmetrizer", "batch_symmetrizer", 
    "remove_equivalent", 
    "vasp_relaxation_settings", "vasp_static_settings", 
    "vasp_launcher", "vasp_batch_launch", 
    "batch_extract_vasp_data", "chgcar_density_switch", "delta_sol_inputs_init", 
    "E_O2_FIT", "U_VALUES", "DELTA_E_M", "EXP_DELTA_H", "EL_PER_XC_VOL", 
    "add_new_dir", "batch_add_new_dirs", 
    "get_elements_from_entries", "init_entries_and_group_by_dim_and_comp", 
    "get_sub_entries", "phase_diagram_init", 
    "calculate_instability_energies", "batch_calculate_instability_energies", 
    "calculate_delta_sol_band_gap", "batch_calculate_delta_sol_band_gaps",
    ]

from pathlib import Path
from typing import Union, Sequence, List, Literal
import itertools
from ruamel.yaml import YAML

PathLike = Union[Path, str]

PMGRelaxSet = Literal[
    'MITRelaxSet', 
    'MPRelaxSet', 
    'MPScanRelaxSet', 
    'MPHSERelaxSet', 
    'MPMetalRelaxSet', 
    'MVLRelax52Set', 
    'MVLScanRelaxSet'
]

PMGStaticSet = Literal[
    'MPStaticSet', 
    'MatPESStaticSet', 
    'MPScanStaticSet'
]

def is_float(string: str) -> bool:
    str_list = string.split(sep='.')
    return ('.' in string) and (len(str_list) <= 2) and all([nbr.isdecimal() for nbr in str_list])


def flatten(sequence: Sequence, level_of_flattening: int = 1) -> List:
    '''
    Unpacks a nested sequence without modifying elements order.

    Parameters:
        sequence (Sequence):        The iterable to unpack.

        level_of_flattening (Int):  The number of nested levels to unpack.
                                    Defaults to 1.
    
    Returns:
        A list flattened the specified number of times.
    '''

    for _ in range(1, level_of_flattening + 1):
        sequence = list(itertools.chain.from_iterable(sequence))
    return sequence


def _yaml_loader(file_path: PathLike):
    if not isinstance(file_path, (Path, str)):
        raise TypeError(f"Expected a Path object or str, got {type(file_path)} instead.")
    
    assert str(file_path).endswith('.yaml'), \
    f'{str(file_path)} is not a .yaml file format.'

    file_path = Path(file_path)

    assert file_path.is_file(), \
    f'{str(file_path)}: no such file found.'

    yaml = YAML()
    with open(file_path, encoding="utf-8") as yaml_file:
        try: yaml_data = yaml.load(yaml_file) or {}
        except Exception as exc:
            warn_msg = f'An exception was thrown during yaml loading of file {str(file_path)}.\n\
                        Data written in this file is ignored to proceed.\n\
                        Thrown exception below:\n\
                        {exc}'
            print(warn_msg)
            return {}
        return dict(yaml_data)

