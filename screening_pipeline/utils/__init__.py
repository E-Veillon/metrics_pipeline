"""Screening pipeline utils module init."""

# I/O modules
from .cif_io import read_cif, write_cif
from .file_io import add_new_dir, batch_add_new_dirs, yaml_loader
from .vasp_io import (
    vasp_relaxation_settings,
    vasp_static_settings,
    vasp_launcher,
    vasp_batch_launch,
    converged_vasprun,
    batch_extract_vasp_data,
    batch_extract_vasp_structures,
    delta_sol_inputs_init,
    extract_vasp_data_for_delta_sol_init,
    dsol_calc_init,
)

# Chemistry related modules
from .fitted_values import E_O2_FIT, U_VALUES, DELTA_E_M, EXP_DELTA_H, EL_PER_XC_VOL
from .spacegroup import structure_symmetrizer, batch_symmetrizer
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

from .matcher import check_interatomic_distances, remove_equivalent, group_by_composition


from .convex_hulls import (
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
    DeltaSolStaticSet,
    get_dsol_struct_dir,
    calc_idx_to_dir_name,
    get_dsol_n_ratio,
    _match_calc_index,
    dsol_calc_init,
    get_dsol_band_gap,
    batch_get_dsol_band_gaps,
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
from .flattener import flatten

__all__ = [
    "has_rare_gas", "discard_rare_gas_structures",
    "has_rare_earth", "discard_rare_earth_structures",
    "get_elements", "get_elemental_subsets", "get_all_elements_groups",
    "get_element_valence_electrons", "get_all_valence_electrons",
    "read_cif", "write_cif",
    "structure_symmetrizer",
    "batch_symmetrizer",
    "remove_equivalent",
    "vasp_relaxation_settings",
    "vasp_static_settings",
    "vasp_launcher",
    "vasp_batch_launch",
    "batch_extract_vasp_data",
    "delta_sol_inputs_init",
    "E_O2_FIT",
    "U_VALUES",
    "DELTA_E_M",
    "EXP_DELTA_H",
    "EL_PER_XC_VOL",
    "add_new_dir",
    "batch_add_new_dirs",
    "get_elements_from_entries",
    "get_sub_entries",
    "get_dsol_band_gap",
    "batch_get_dsol_band_gaps",
    "load_phase_diagram_entries",
    "yaml_loader",
]
