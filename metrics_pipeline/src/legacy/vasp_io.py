"""Legacy functions from vasp_io.py module in versions < 2.0.0, to be replaced later."""

import os
import json
import typing as tp
import warnings

from monty.dev import deprecated
import xml.etree.ElementTree as ET

from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun

from metrics_pipeline.src.io import PathLike, check_file_or_dir


def _converged_vasprun(run_path: PathLike, **kwargs) -> tp.Union[Vasprun, None]:
    """
    Read and parse vasprun.xml file at given VASP run directory.
    Check whether the VASP run terminated and converged normally.
    Returns a Vasprun object if it is the case, or None otherwise.
    
    Parameters:
        run_path (Path|str):    Path to the directory containing the VASP calculation.

        **kwargs:               Keyword arguments supported by pymatgen Vasprun class.
    
    Returns:
        Vasprun object if the run terminated and converged normally,
        None otherwise.
    """
    check_file_or_dir(run_path, "dir")
    vasprun_path = os.path.join(run_path, "vasprun.xml")
    check_file_or_dir(vasprun_path, "file", allowed_formats="xml")

    try:
        vasprun = Vasprun(filename=vasprun_path, **kwargs)
    except ET.ParseError:
        return None
    except UnicodeDecodeError:
        warnings.warn(
            f"WARNING: vasprun.xml file at {run_path} contains "
            "unreadable characters for 'utf-8' codec.\n"
            "Associated data is therefore considered erroneous "
            "and is not parsed further."
        )
        return None

    if not vasprun.converged:
        return None

    return vasprun


def _check_summary_data(
    struct_dir: PathLike, summary_file: PathLike, key_to_check: str
) -> bool:
    """
    Check whether given structure directory is present in the summary file and
    whether it was rejected in the previous step.
    """
    struct_dir = str(struct_dir)

    with open(summary_file, "rt", encoding="utf-8") as fp:
        summary = json.load(fp)

    try:
        prev_struct_data = next(filter(
            lambda data: os.path.samefile(data["path"], struct_dir),
            summary
        ))
    except StopIteration:
        warnings.warn(
            f"Structure path '{struct_dir}' was not found in the summary file.\n"
            "You may want to pass it to the previous screening step before this one.\n"
            "This structure is assumed not viable and is ignored for this step.",
            stacklevel=2
        )
        return False

    return prev_struct_data[key_to_check]


def _check_vasp_data(
    struct_dir: PathLike,
    path_to_summary: tp.Optional[PathLike] = None,
    key_to_check: tp.Optional[str] = None
) -> Vasprun|None:
    """
    Check whether given path leads to a previously accepted structure,
    its vasprun.xml file exists and is a normally terminated and converged run.

    Parameters:
        struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.
        
        path_to_summary (str):  Path to a JSON summary file containing results from previous step.
                                Used to not consider structures that failed before.
                                If not provided, all structure subdirs will be extracted.

        key_to_check (str):     The dict key associated to the bool used to verify eligibility.
                                If path_to_summary is given, it must be given too.

    Returns: Vasprun object if it is eligible and normally converged, None otherwise.
    """
    check_file_or_dir(struct_dir, "dir")
    struct_dir = str(struct_dir)

    if path_to_summary is not None and key_to_check is None:
        raise ValueError("If path_to_summary is given, key_to_check must be given too.")

    if (
        path_to_summary is not None
        and not _check_summary_data(struct_dir, path_to_summary, key_to_check) # type: ignore
    ):
        return None

    return _converged_vasprun(
        struct_dir, parse_dos=False, parse_eigen=False, parse_potcar_file=False
    )


def extract_vasp_data_for_delta_sol_init(
    struct_dir: PathLike,
    path_to_summary: tp.Optional[PathLike] = None,
    key_to_check: tp.Optional[str] = None
) -> tuple[str, dict[str, tp.Union[Structure, float]]]:
    """
    Extracts VASP data from a previous run for one structure.
    Keeps only relevant data for Δ-Sol method.

    Parameters:
        struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.
        
        path_to_summary (str):  Path to a JSON summary file containing results from previous step.
                                Used to not consider structures that failed before.
                                If not provided, all structure subdirs will be extracted.

        key_to_check (str):     The dict key associated to the bool used to verify eligibility.
                                If path_to_summary is given, it must be given too.

    Returns:
        Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict,
                                containing following data, used for Δ-Sol method:
                                    - structure itself,
                                    - its final energy (in eV).
    """
    vasprun = _check_vasp_data(struct_dir, path_to_summary, key_to_check)

    if vasprun is None:
        return "", {}

    struct_name = os.path.basename(str(struct_dir))
    structure = vasprun.final_structure
    final_energy = vasprun.final_energy

    struct_dict = {
        "structure": structure,
        "final_energy": final_energy,
    }
    struct_data = (struct_name, struct_dict)

    return struct_data