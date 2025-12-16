"""Screening pipeline utils subpackage."""

# I/O package
from .io import (
    PathLike, MAINDIRPATH, SCRIPTSPATH, UTILSPATH, CONFIGPATH,
    check_file_format, check_file_or_dir,
    read_cif, symmetrize_and_write_cif,
    PDDataset, MPDatasetDownloader,
    JsonLoader, JsonWriter,
    PoscarBlock, PoscarFile,
    VaspWriter, VaspParser, VaspExtractor, ExtractMethod,
    load_yaml_as_dict,
)
from .computations.vasp import (
    PMGRelaxSet, PMGStaticSet, vasp_relaxation_settings, vasp_static_settings, dsol_calc_init
)

# Metrics package
from .metrics import (
    StructValidity, Viability, Symmetry, ElementaryMetastability,
    Stability, Unicity, Novelty, SUN,
    Coverage, EMD, FrechetDistance, RMSD
)

# Chemistry related modules
from .utils import (
    FormulaLike,
    has_rare_gas, discard_rare_gas_structures,
    has_rare_earth, discard_rare_earth_structures,
    get_elements, get_elemental_subsets,
    get_all_elements_groups,
    get_element_valence_electrons, get_all_valence_electrons,
)
from .legacy.fitted_values import U_VALUES, EL_PER_XC_VOL
from .legacy.spacegroup import structure_symmetrizer, batch_symmetrizer, SPG_NUM_TO_PG, PG_TO_SYSTEM

# Pipeline steps modules
from .legacy.delta_sol import (
    DSolStaticSet,
    get_dsol_struct_dir,
    calc_idx_to_dir_name,
    get_dsol_n_ratio,
    get_dsol_band_gap,
    batch_get_dsol_band_gaps,
)

# Metrics related modules
from .computations.models import vectors_from_alignn, get_crystalnn_fingerprints
from .computations.local import get_densities
from .legacy.matcher import (
    check_viability,
    check_interatomic_distances, group_by_composition,
    batch_group_by_equivalence, remove_equivalent,
    batch_get_novel_structures,
    is_struct_dir, match_struct_dirs
)

# Other modules
from .utils import (
    check_type, check_num_value, flatten, VisualIterator, redirect_c_stdout, redirect_c_stderr
)

__all__ = [
    # I/O package
    "PathLike", "MAINDIRPATH", "SCRIPTSPATH", "UTILSPATH", "CONFIGPATH",
    "check_file_format", "check_file_or_dir",
    "read_cif", "symmetrize_and_write_cif",
    "PoscarBlock", "PoscarFile",
    "VaspWriter", "VaspParser", "VaspExtractor", "ExtractMethod",
    "PMGRelaxSet", "PMGStaticSet",
    "vasp_relaxation_settings", "vasp_static_settings",
    "dsol_calc_init",
    "load_yaml_as_dict",
    # Metrics package
    "StructValidity", "Viability", "Symmetry", "ElementaryMetastability",
    "Stability", "Unicity", "Novelty", "SUN",
    "Coverage", "EMD", "FrechetDistance", "RMSD",
    # Chemistry
    "U_VALUES", "EL_PER_XC_VOL",    # fitted_values
    "structure_symmetrizer", "batch_symmetrizer",   # spacegroup
    "FormulaLike", "has_rare_gas", "discard_rare_gas_structures",       # periodic_table
    "has_rare_earth", "discard_rare_earth_structures",                  #
    "get_elements", "get_elemental_subsets", "get_all_elements_groups", #
    "get_element_valence_electrons", "get_all_valence_electrons",       #
    # Pipeline steps
    "DSolStaticSet",            # delta_sol
    "get_dsol_struct_dir",      #
    "calc_idx_to_dir_name",     #
    "get_dsol_n_ratio",         #
    "dsol_calc_init",           #
    "get_dsol_band_gap",        #
    "batch_get_dsol_band_gaps", #
    # Metrics
    "vectors_from_alignn",
    "get_densities",    # density
    "check_viability", "check_interatomic_distances", "group_by_composition",   # matcher
    "batch_group_by_equivalence", "remove_equivalent",                          #
    "batch_get_novel_structures",                                               #
    "is_struct_dir", "match_struct_dirs",                                       #
    # Other
    "check_type", "check_num_value",    # common_asserts
    "flatten",  # flattener
    "redirect_c_stdout", "redirect_c_stderr",   # redirect
    "VisualIterator",   # visual_iterator
]
