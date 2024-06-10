from .periodic_table import (
    has_rare_gas,
    discard_rare_gas_structures,
    has_rare_earth,
    discard_rare_earth_structures,
    get_elements,
    get_elemental_subsets,
    get_all_elements_groups,
    get_element_valence_electrons,
    get_all_valence_electrons,
)
from .cif_io import read_cif, write_cif, extract_cif_from_file, cif_str_to_struct
from .matcher import (
    flatten,
    remove_equivalent,
    group_by_composition,
    batch_group_by_equivalence,
    batch_get_novel_structures
)
from .vasp_io import (
    vasp_relaxation_settings,
    vasp_static_settings,
    vasp_launcher,
    vasp_batch_launch,
    converged_Vasprun,
    batch_extract_vasp_data,
    batch_extract_vasp_structures,
    delta_sol_inputs_init,
    extract_vasp_data_for_delta_sol_init,
    delta_sol_calculation_init,
)
from .spacegroup import structure_symmetrizer, batch_symmetrizer
from .fitted_values import E_O2_FIT, U_VALUES, DELTA_E_M, EXP_DELTA_H, EL_PER_XC_VOL
from .paths import add_new_dir, batch_add_new_dirs
from .data_process import (
    check_interatomic_distances,
    get_max_dim,
    init_entries_from_dict,
    filter_database_entries,
    get_elements_from_entries,
    group_by_dim_and_comp,
    #init_entries_and_group_by_dim_and_comp,
    get_sub_entries,
    get_lacking_elts_entries,
    batch_compute_e_above_hull,
    calculate_delta_sol_band_gap,
    batch_calculate_delta_sol_band_gaps,
)
from .custom_types import (
    PathLike, FormulaLike, PMGRelaxSetType, PMGStaticSetType,
    PMGRelaxSet, PMGStaticSet
)
from .dataset import load_phase_diagram_entries, mp_api_download, process_oqmd_json_file
from .ml_vectors import vectors_from_alignn
from .distribution import (
        recall,
        precision,
        frechet_distance,
        wasserstein_distance,
)
from .density import get_densities
from .rmsd import rmsd_from_structures
from .crystalnn import to_crystalnn_fingerprint
from .utils import _yaml_loader

__all__ = [
    "has_rare_gas",
    "discard_rare_gas_structures",
    "has_rare_earth",
    "discard_rare_earth_structures",
    "get_elements",
    "get_elemental_subsets",
    "get_all_elements_groups",
    "get_all_valence_electrons",
    "read_cif",
    "write_cif",
    "structure_symmetrizer",
    "batch_symmetrizer",
    "remove_equivalent",
    "vasp_relaxation_settings",
    "vasp_static_settings",
    "vasp_launcher",
    "vasp_batch_launch",
    "batch_extract_vasp_data",
    "chgcar_density_switch",
    "delta_sol_inputs_init",
    "E_O2_FIT",
    "U_VALUES",
    "DELTA_E_M",
    "EXP_DELTA_H",
    "EL_PER_XC_VOL",
    "add_new_dir",
    "batch_add_new_dirs",
    "get_elements_from_entries",
    "init_entries_and_group_by_dim_and_comp",
    "get_sub_entries",
    "phase_diagram_init",
    "calculate_instability_energies",
    "batch_calculate_instability_energies",
    "calculate_delta_sol_band_gap",
    "batch_calculate_delta_sol_band_gaps",
    "load_phase_diagram_entries",
]
