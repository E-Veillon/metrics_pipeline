"""Package containing all I/O operations."""
# TODO: split some functionalities to actually have a pure I/O package.
from .io_base import (
    PathLike, MAINDIRPATH, SCRIPTSPATH, UTILSPATH, CONFIGPATH, check_file_format, check_file_or_dir
)
from .cif import read_cif, symmetrize_and_write_cif
from .poscar import PoscarBlock, PoscarFile
from .vasp_files import (
    PMGRelaxSet, PMGStaticSet, write_and_run_vasp,
    vasp_relaxation_settings, vasp_static_settings,
    get_struct_from_vasp, converged_vasprun,
    extract_vasp_data_for_delta_sol_init, batch_extract_vasp_data,
    batch_extract_vasp_structures, dsol_calc_init
)
from .yaml import load_yaml_as_dict


__all__ = [
    "PathLike", "MAINDIRPATH", "SCRIPTSPATH", "UTILSPATH", "CONFIGPATH",
    "check_file_format", "check_file_or_dir",
    "read_cif", "symmetrize_and_write_cif",
    "PoscarBlock", "PoscarFile",
    "write_and_run_vasp",
    "vasp_relaxation_settings", "vasp_static_settings",
    "get_struct_from_vasp", "converged_vasprun",
    "extract_vasp_data_for_delta_sol_init", "batch_extract_vasp_data",
    "batch_extract_vasp_structures", "dsol_calc_init",
    "load_yaml_as_dict",
]