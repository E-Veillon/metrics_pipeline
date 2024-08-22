"""Screening pipeline utils module init."""

# I/O modules
from .cif_io import read_cif, write_cif
from .file_io import (
    MAINDIRPATH, SCRIPTSPATH, UTILSPATH, CONFIGPATH,
    check_file_format, check_file_or_dir, add_new_dir, yaml_loader

)
from .vasp_io import (
    write_and_run_vasp,
    vasp_relaxation_settings, vasp_static_settings,
    get_struct_from_vasp, converged_vasprun,
    extract_vasp_data_for_delta_sol_init, batch_extract_vasp_data,
    batch_extract_vasp_structures, dsol_calc_init
)

# Chemistry related modules
from .fitted_values import U_VALUES, EL_PER_XC_VOL
from .spacegroup import structure_symmetrizer, batch_symmetrizer
from .periodic_table import (
    has_rare_gas, discard_rare_gas_structures,
    has_rare_earth, discard_rare_earth_structures,
    get_elements, get_elemental_subsets,
    get_all_elements_groups,
    get_element_valence_electrons, get_all_valence_electrons,
)

# Pipeline steps modules
from .phase_diagrams import (
    get_max_dim,
    init_entries_from_dict,
    filter_database_entries,
    get_elements_from_entries,
    group_by_dim_and_comp,
    get_sub_entries,
    get_lacking_elts_entries,
    batch_compute_e_above_hull,
)
from .delta_sol import (
    DSolStaticSet,
    get_dsol_struct_dir,
    calc_idx_to_dir_name,
    get_dsol_n_ratio,
    get_dsol_band_gap,
    batch_get_dsol_band_gaps,
)

# Metrics related modules
from .crystalnn import to_crystalnn_fingerprint
from .dataset import (
    load_phase_diagram_entries,
    mp_api_download, process_oqmd_json_file
)
from .density import get_densities
from .distribution import (
    recall, precision, frechet_distance, get_emd
)
from .matcher import (
    check_interatomic_distances, group_by_composition,
    batch_group_by_equivalence, remove_equivalent,
    batch_get_novel_structures,
    is_struct_dir, match_struct_dirs
)
from .ml_vectors import vectors_from_alignn
from .rmsd import rmsd_from_structures

# Other modules
from .common_asserts import check_type, check_num_value
from .custom_types import (
    PathLike, FormulaLike,
    PMGRelaxSet, PMGStaticSet
)
from .flattener import flatten
from .redirect import redirect_c_stdout, redirect_c_stderr

__all__ = [
    # I/O
    "read_cif", "write_cif",    # cif_io
    "MAINDIRPATH", "SCRIPTSPATH", "UTILSPATH", "CONFIGPATH",                # file_io
    "check_file_format", "check_file_or_dir", "add_new_dir", "yaml_loader", #
    "write_and_run_vasp",
    "vasp_relaxation_settings", "vasp_static_settings",                 # vasp_io
    "get_struct_from_vasp", "converged_vasprun",                        #
    "extract_vasp_data_for_delta_sol_init", "batch_extract_vasp_data",  #
    "batch_extract_vasp_structures", "dsol_calc_init",                  #
    # Chemistry
    "U_VALUES", "EL_PER_XC_VOL",    # fitted_values
    "structure_symmetrizer", "batch_symmetrizer",   # spacegroup
    "has_rare_gas", "discard_rare_gas_structures",                      # periodic_table
    "has_rare_earth", "discard_rare_earth_structures",                  #
    "get_elements", "get_elemental_subsets", "get_all_elements_groups", #
    "get_element_valence_electrons", "get_all_valence_electrons",       #
    # Pipeline steps
    "get_max_dim",                  #convex_hulls
    "init_entries_from_dict",       #
    "filter_database_entries",      #
    "get_elements_from_entries",    #
    "group_by_dim_and_comp",        #
    "get_sub_entries",              #
    "get_lacking_elts_entries",     #
    "batch_compute_e_above_hull",   #
    "DSolStaticSet",            # delta_sol
    "get_dsol_struct_dir",      #
    "calc_idx_to_dir_name",     #
    "get_dsol_n_ratio",         #
    "dsol_calc_init",           #
    "get_dsol_band_gap",        #
    "batch_get_dsol_band_gaps", #
    # Metrics
    "to_crystalnn_fingerprint", # crystalnn
    "load_phase_diagram_entries",                   # dataset
    "mp_api_download", "process_oqmd_json_file",    #
    "get_densities",    # density
    "recall", "precision", "frechet_distance", "get_emd",  # distribution
    "check_interatomic_distances", "group_by_composition",  # matcher
    "batch_group_by_equivalence", "remove_equivalent",      #
    "batch_get_novel_structures",                           #
    "is_struct_dir", "match_struct_dirs",                   #
    "vectors_from_alignn",  # ml_vectors
    "rmsd_from_structures", # rmsd
    # Other
    "check_type", "check_num_value",    # common_asserts
    "PathLike", "FormulaLike",              # custom_types
    "PMGRelaxSet", "PMGStaticSet",          #
    "flatten",  # flattener
    "redirect_c_stdout", "redirect_c_stderr",   # redirect
]
