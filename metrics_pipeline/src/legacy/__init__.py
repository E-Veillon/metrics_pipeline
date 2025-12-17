"""
Subpackage containing legacy modules from versions <2.0.0 that are kept temporarily
until a better replacement is made.
"""
import warnings
warnings.warn(
    "Use of the legacy package is deprecated. "
    "Removal of its features may happen as soon as an update implements better replacements.",
    DeprecationWarning
)

from .delta_sol import (
    EL_PER_XC_VOL, DSolCalc, get_dsol_struct_dir, calc_idx_to_dir_name,
    get_dsol_n_ratio, get_dsol_band_gap, batch_get_dsol_band_gaps, 
)
from .fitted_values import E_O2_FIT, EXP_DELTA_H, DELTA_E_M
from .matcher import batch_group_by_equivalence, remove_equivalent
from .vasp_io import extract_vasp_data_for_delta_sol_init