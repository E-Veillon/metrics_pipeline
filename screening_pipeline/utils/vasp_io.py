#!/usr/bin/python
"""
Functions to manage, read, and write VASP files and launch VASP computations.
"""


import os
import json
import warnings
from functools import partial
from typing import Optional, Dict, List, Union, Sequence, Tuple, Literal, Any
from pathlib import Path
import xml.etree.ElementTree as ET
from monty.os import cd
from tqdm.contrib.concurrent import process_map

# PYTHON MATERIAL GENOMICS
from pymatgen.core import Structure, SiteCollection
from pymatgen.io.vasp import VaspInput, Vasprun, Poscar, Xdatcar
from pymatgen.io.vasp.sets import (
    DictSet, MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPHSERelaxSet,
    MPMetalRelaxSet, MVLRelax52Set, MVLScanRelaxSet,
    MPStaticSet, MatPESStaticSet, MPScanStaticSet, MPSOCSet
)

# LOCAL IMPORTS
from .common_asserts import check_type, check_num_value
from .custom_types import PathLike, PMGRelaxSet, PMGStaticSet
from .delta_sol import DSolStaticSet, get_all_valence_electrons, get_dsol_n_ratio
from .fitted_values import U_VALUES
from .file_io import check_file_or_dir
from .matcher import match_struct_dirs


########################################


def _mitrelaxset_incar_corrections(n_sites: int|None = None) -> Dict:
    """
    Systematic correction for MITRelaxSet INCAR tags that do not match with 
    parameters given in the original work of Jain et al.

    Reference:
        A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller,
        K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).

    Parameters:
        n_sites (int):  The number of sites in the structure. Used to pass explicitly
                        an EDIFF tag to pymatgen. Usually it is not necessary, as
                        the EDIFF_PER_ATOM tag is supported for the same function.
                        If not given, the correction will replace the default EDIFF
                        by an EDIFF_PER_ATOM of 5e-5 eV/atom.

    Returns:
        A dictionnary containing the INCAR tags corrections for MITRelaxSet.
    """
    corrections = {}

    if n_sites is None:
        corrections.update(
            {
                "EDIFF": None,
                "EDIFF_PER_ATOM": round(float(5e-5), 6)
            }
        )
    else:
        corrections.update({"EDIFF": round(float(5e-5) * n_sites, 6)})

    corrections.update(
        {
            "LDAUU": U_VALUES,
            "LDAUL": {
                "F": {
                    "Ag": 2, "Co": 2, "Cr": 2, "Cu": 2, "Fe": 2,
                    "Mn": 2, "Mo": 2, "Nb": 2, "Ni": 2, "Re": 2,
                    "Ta": 2, "V": 2, "W": 2
                },
                "O": {
                    "Ag": 2, "Co": 2, "Cr": 2, "Cu": 2, "Fe": 2,
                    "Mn": 2, "Mo": 2, "Nb": 2, "Ni": 2, "Re": 2,
                    "Ta": 2, "V": 2, "W": 2
                },
                "S": {
                    "Fe": 2,
                    "Mn": 2,  #"Mn": 2.5 -> 2 (quantum number l have to be an integer)
                },
            }
        }
    )
    return corrections


########################################


def _check_vasp_input(vasp_input: VaspInput) -> None:
    """Perform several tests to verify VaspInput correctness."""
    check_type(vasp_input, "vasp_input", (VaspInput,))
    if vasp_input.get("INCAR") is None:
        raise ValueError(
            "check_vasp_input: There is no INCAR defined in the input !"
        )
    if vasp_input.get("POSCAR") is None:
        raise ValueError(
            "check_vasp_input: There is no POSCAR defined in the input !"
        )
    if (
        vasp_input.get("KPOINTS") is None
        and vasp_input["INCAR"].get("KSPACING") is None
    ):
        raise ValueError(
            "check_vasp_input: There is no KPOINTS or KSPACING tag defined in the input !"
        )
    if vasp_input.get("POTCAR") is None:
        raise ValueError(
            "check_vasp_input: There is no POTCAR defined in the input !"
        )


########################################


def write_and_run_vasp(
    vasp_input: VaspInput, run_path: PathLike, vasp_exe: PathLike = "vasp"
) -> None:
    """
    Function to run VASP from a VaspInput object.

    Parameters:
        vasp_input (VaspInput): The VaspInput object containing all necessary
                                data to write VASP input files.

        run_path (str|Path):    Path to the directory where VASP files
                                will be written and run.

        vasp_exe (str|Path):    Absolute path to the VASP executable.
                                If not given, attempt the usual shortcut
                                "vasp" launch command at given path.
    """
    _check_vasp_input(vasp_input)
    check_file_or_dir(run_path, "dir")
    if vasp_exe != "vasp":
        check_file_or_dir(vasp_exe, "file")

    vasp_input.write_input(output_dir=run_path)

    with cd(run_path):
        os.system(f"{vasp_exe}")


########################################


def _relax_set_init(
    structure: SiteCollection,
    preset: PMGRelaxSet = PMGRelaxSet.MPRELAXSET,
    corrections: Optional[Dict] = None,
) -> DictSet:
    """Init a relaxation set of VASP input files."""
    corrections = corrections or {}
    incar_corrections = {}

    if preset == PMGRelaxSet.MITRELAXSET:
        incar_corrections = _mitrelaxset_incar_corrections(structure.num_sites)

    incar_corrections.update(corrections.get("INCAR", {}))
    kpoints_corrections = corrections.get("KPOINTS", {})
    potcar_corrections = corrections.get("POTCAR", {})
    potcar_functional_correction = corrections.get("POTCAR_FUNCTIONAL", {})

    match preset.value:
        case "MITRelaxSet":
            preset_obj = MITRelaxSet(
                structure=structure,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MPRelaxSet":
            preset_obj = MPRelaxSet(
                structure=structure,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MPScanRelaxSet":
            preset_obj = MPScanRelaxSet(
                structure=structure,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MPHSERelaxSet":
            preset_obj = MPHSERelaxSet(
                structure=structure,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MPMetalRelaxSet":
            preset_obj = MPMetalRelaxSet(
                structure=structure,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MVLRelax52Set":
            preset_obj = MVLRelax52Set(
                structure=structure,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MVLScanRelaxSet":
            preset_obj = MVLScanRelaxSet(
                structure=structure,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )

    return preset_obj


########################################


def _static_set_init(
    structure: SiteCollection,
    preset: Union[PMGStaticSet, Literal["DsolStaticSet"]] = PMGStaticSet.MPSTATICSET,
    nelect: float|None = None,
    corrections: Optional[Dict] = None,
) -> DictSet:
    """Init a static calculation set of VASP input files."""

    corrections = corrections or {}
    incar_corrections = corrections.get("INCAR", {})
    kpoints_corrections = corrections.get("KPOINTS", {})
    potcar_corrections = corrections.get("POTCAR", {})
    potcar_functional_correction = corrections.get("POTCAR_FUNCTIONAL", {})

    if preset == "DSolStaticSet":
        preset_obj = DSolStaticSet(
            structure=structure,
            incar_nelect=nelect,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    else:
        match preset.value:
            case "MPStaticSet":
                preset_obj = MPStaticSet(
                    structure=structure,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )
            case "MatPESStaticSet":
                preset_obj = MatPESStaticSet(
                    structure=structure,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )
            case "MPScanStaticSet":
                preset_obj = MPScanStaticSet(
                    structure=structure,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )
            case "MPSOCSet":
                preset_obj = MPSOCSet(
                    structure=structure,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )

    return preset_obj


########################################


def vasp_relaxation_settings(
    structure: SiteCollection,
    preset: PMGRelaxSet = PMGRelaxSet.MPRELAXSET,
    user_corrections: Optional[Dict] = None,
) -> VaspInput:
    """
    Setup VASP inputs for a given structure using one of the pymatgen relaxation presets.

    Parameters:
        structure (SiteCollection):     The structure to write VASP inputs for.

        preset (str):                   The pymatgen preset to use for VASP inputs
                                        initialization.

        user_corrections (dict):        User defined settings. It allows to override
                                        some of the preset INCAR, KPOINTS or POTCAR
                                        settings if necessary. Defaults to None.
    """
    check_type(structure, "structure", (SiteCollection,))
    assert preset in PMGRelaxSet, (
        "'preset' argument not recognized. "
        "It must be one of the allowed pymatgen relaxation presets:\n"
        f"{PMGRelaxSet.values}"
    )
    if user_corrections is not None:
        check_type(user_corrections, "user_corrections", (Dict,))

    vasp_input = _relax_set_init(
        structure, preset, user_corrections
    ).get_input_set()

    return vasp_input


########################################


def vasp_static_settings(
    structure: Optional[SiteCollection] = None,
    preset: Union[PMGStaticSet, Literal["DsolStaticSet"]] = PMGStaticSet.MPSTATICSET,
    nelect: float|None = None,
    user_corrections: Optional[dict] = None,
) -> VaspInput:
    """
    Setup VASP inputs for a given structure using one of the pymatgen static presets.

    Parameters:
        structure (SiteCollection):     The structure to write VASP inputs for.

        preset (str):                   The pymatgen preset to use for VASP inputs initialization.
                                        Can also be the homemade "DSolStaticSet" if Δ-Sol
                                        method by Chan et al. (2010) is used.

        nelect (float):                 Only useful if DSolStaticSet is used.
                                        Sets the NELECT tag in INCAR file.
                                        Ignored if from_prev_calc is True.

        user_corrections (dict):        User defined settings. It allows to override some of 
                                        the preset INCAR, KPOINTS or POTCAR settings if necessary.
                                        Defaults to None.
    """
    if structure is not None:
        check_type(structure, "structure", (SiteCollection,))
    assert preset in PMGStaticSet, (
        "'preset' argument not recognized. "
        "It must be one of the allowed pymatgen static presets:\n"
        f"{PMGStaticSet.values}."
    )
    if user_corrections is not None:
        check_type(user_corrections, "user_corrections", (Dict,))

    vasp_input = _static_set_init(
        structure, preset, nelect, user_corrections
    ).get_input_set()

    return vasp_input


########################################


def get_struct_from_vasp(
    path: PathLike,
    try_vasprun: bool = True,
    try_contcar: bool = True,
    try_xdatcar: bool = True,
    try_poscar: bool = True
) -> Structure:
    """
    Attempt to extract the best structure possible from a VASP run directory.
    
    Algorithm:
        Attempt to extract a structure from a VASP run in following order,
        passing to the next step if the current one fails:
        1 - Tries to get the final structure from vasprun.xml
        2 - Tries to get the final structure from CONTCAR
        3 - Tries to get the last structure from XDATCAR
        4 - Tries to get the initial structure from POSCAR
    
    Parameters:
        path (str|Path): Path to the VASP run directory.

        try_vasprun (bool): Whether to try to get a structure from vasprun.xml.
                            Defaults to True.

        try_contcar (bool): Whether to try to get a structure from CONTCAR.
                            Defaults to True.

        try_xdatcar (bool): Whether to try to get a structure from XDATCAR.
                            Defaults to True.

        try_poscar (bool):  Whether to try to get a structure from POSCAR.
                            Defaults to True.

    Raises:
        FileNotFoundError if none of the attempts were successful.

    Returns:
        A Structure object if it could be extracted.
    """
    check_file_or_dir(path, "dir")

    if try_vasprun:
        try:
            file = os.path.join(path, "vasprun.xml")
            structure = Vasprun(
                file, parse_dos=False, parse_eigen=False, parse_potcar_file=False
            ).final_structure
            return structure
        except (FileNotFoundError, ET.ParseError, UnicodeDecodeError):
            pass

    if try_contcar:
        try:
            file = os.path.join(path, "CONTCAR")
            structure = Poscar.from_file(file).structure
            return structure
        except FileNotFoundError:
            pass

    if try_xdatcar:
        try:
            file = os.path.join(path, "XDATCAR")
            structure = Xdatcar(file).structures[-1]
            return structure
        except FileNotFoundError:
            pass

    if try_poscar:
        try:
            file = os.path.join(path, "POSCAR")
            structure = Poscar.from_file(file).structure
            return structure
        except FileNotFoundError:
            pass

    raise FileNotFoundError(
        f"get_struct_from_vasp: Unable to extract a structure from '{path}'.\n"
        "Tried files:\n"
        f"- vasprun.xml: {try_vasprun}\n"
        f"- CONTCAR: {try_contcar}\n"
        f"- XDATCAR: {try_xdatcar}\n"
        f"- POSCAR: {try_poscar}\n"
    )


########################################


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
            stack_level=2
        )
        return False

    return prev_struct_data[key_to_check]


########################################

def converged_vasprun(run_path: PathLike, **kwargs) -> Union[Vasprun, None]:
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


########################################


def _check_vasp_data(
    struct_dir: PathLike,
    path_to_summary: Optional[PathLike] = None,
    key_to_check: Optional[str] = None
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
    struct_dir: str = str(struct_dir)

    if (
        path_to_summary is not None
        and not _check_summary_data(struct_dir, path_to_summary, key_to_check)
    ):
        return None

    return converged_vasprun(
        struct_dir, parse_dos=False, parse_eigen=False, parse_potcar_file=False
    )


########################################


def extract_vasp_data_for_convex_hull(
    struct_dir: PathLike,
    path_to_summary: Optional[PathLike] = None,
    key_to_check: Optional[str] = None
) -> Tuple[str, Dict[str, Union[Structure, float]]]:
    """
    Extracts VASP data from a previous run for one structure.
    Keeps only relevant data for relative stability calculation.

    Parameters:
        struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.
        
        path_to_summary (str):  Path to a JSON summary file containing results from previous step.
                                Used to not consider structures that failed before.
                                If not provided, all structure subdirs will be extracted.

        key_to_check (str):     The dict key associated to the bool used to verify eligibility.
                                If path_to_summary is given, it must be given too.

    Returns:
        Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict,
                                containing following data, used in stability calculation:
                                    - structure directory name,
                                    - structure chemical composition as Composition object,
                                    - generated structure energy (in eV),
    """
    vasprun = _check_vasp_data(struct_dir, path_to_summary, key_to_check)

    if vasprun is None:
        return {}

    # We want data of the generated structure for convex hulls, not the relaxed one
    struct_name  = os.path.basename(str(struct_dir))
    composition = vasprun.initial_structure.composition
    generated_energy: float = vasprun.ionic_steps[0]["e_0_energy"]


    struct_dict = {
        "entry_id": struct_name,
        "composition": composition,
        "final_energy": generated_energy
    }
    struct_data = (struct_name, struct_dict)

    return struct_data


########################################


def extract_vasp_data_for_delta_sol_init(
    struct_dir: PathLike,
    path_to_summary: Optional[PathLike] = None,
    key_to_check: Optional[str] = None
) -> Tuple[str, Dict[str, Union[Structure, float]]]:
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
        return {}

    struct_name = os.path.basename(str(struct_dir))
    structure = vasprun.final_structure
    final_energy = vasprun.final_energy

    struct_dict = {
        "structure": structure,
        "final_energy": final_energy,
    }
    struct_data = (struct_name, struct_dict)

    return struct_data


########################################


def extract_vasp_data_for_delta_sol_calc(
    struct_dir: PathLike,
    path_to_summary: Optional[PathLike] = None,
    key_to_check: Optional[str] = None
) -> Tuple[str, Dict[str, Union[Structure, float]]]:
    """
    Extract the results of Δ-Sol computations.
    
    Parameters:
        struct_dir (str|Path):  Structure directory containing subdirs of Δ-Sol computations.
        
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
    check_file_or_dir(struct_dir, "dir")

    calc_dirs   = list(filter(os.path.isdir, os.listdir(struct_dir)))
    struct_dict = {}

    for calc_dir in calc_dirs:
        if (
            path_to_summary is not None
            and not _check_summary_data(struct_dir, path_to_summary, key_to_check)
        ):
            return {}

        calc_data = extract_vasp_data_for_delta_sol_init(
            calc_dir, path_to_summary, key_to_check
        )
        if not calc_data:
            msg = f"WARNING: {calc_dir} could not be parsed, either because the "
            msg += "VASP run terminated on an error, on a timeout limit "
            msg += "or it did not converge after the maximum ionic step was reached."
            warnings.warn(msg)
            continue

        if struct_dict.get("structure") is None:
            struct_dict.update({"structure": calc_data[1]["structure"]})

        struct_dict.update({calc_data[0]: calc_data[1]["final_energy"]})

    struct_name = Path(struct_dir).name
    struct_data = (struct_name, struct_dict)

    return struct_data


########################################


def batch_extract_vasp_data(
        method: Literal["convex_hull", "delta_sol_init", "delta_sol_calc"],
        base_dir: PathLike,
        structs_names: Optional[Sequence[str]] = None,
        path_to_summary: Optional[PathLike] = None,
        key_to_check: Optional[str] = None,
        workers: int = 1
) -> Dict[str, Dict[str, Any]]:
    """
    Extracts VASP data from a previous run for each structure directory in given directory.
    Keeps only data that are useful according to given method.

    Parameters:
        method (str):           Name of the method that will use the data, used to know
                                which data should be extracted.
                                Actual methods supported: 'convex_hull', 'delta_sol_init',
                                'delta_sol_calc'.

        base_dir (str|Path):    Directory containing structures subdirs to extract data from.

        structs_names ([str]):  Provide specific structures sub-directories to extract data from.
                                If specified, only specified subdirs in base_dir are checked.
                                If not, all subdirs in base_dir are checked.
        
        path_to_summary (str):  Path to a JSON summary file containing results from previous steps.
                                Used to not consider structures that failed before.
                                If not provided, all structure subdirs will be extracted.

        key_to_check (str):     The dict key associated to the bool used to verify eligibility.
                                If path_to_summary is given, it must be given too.

        workers (int):          Number of parallel processes to spawn.

    Returns:
        Dict[str, Dict]:        Dict with structure directory names as keys,
                                and a dict containing useful data according to chosen method
                                for corresponding structure as values.

                                Data returned for 'convex_hull' method:
                                    - composition of the formula unit,
                                    - energy of the unrelaxed structure in eV.

                                Data returned for 'delta_sol' method:
                                    - structure itself,
                                    - final energy of the relaxation in eV.
    """
    check_file_or_dir(base_dir, "dir")
    assert os.listdir(str(base_dir)), (
        f"{str(base_dir)}: Directory exists but is empty."
    )
    if path_to_summary is not None:
        check_type(path_to_summary, "path_to_summary", (Path, str))
    check_type(workers, "workers", (int,))
    check_num_value(workers, "workers", ">", 0)

    match method:
        case "convex_hull":
            set_vasp_extractor = partial(
                extract_vasp_data_for_convex_hull,
                path_to_summary=path_to_summary,
                key_to_check=key_to_check
            )
        case "delta_sol_init":
            set_vasp_extractor = partial(
                extract_vasp_data_for_delta_sol_init,
                path_to_summary=path_to_summary,
                key_to_check=key_to_check
            )
        case "delta_sol_calc":
            set_vasp_extractor = partial(
                extract_vasp_data_for_delta_sol_calc,
                path_to_summary=path_to_summary,
                key_to_check=key_to_check
            )
        case str():
            raise NotImplementedError(
                f"Provided method ({method}) is not supported.\n"
                "Supported methods are: "
                "'convex_hull', 'delta_sol_init', 'delta_sol_calc'."
            )
        case _:
            check_type(method, "method", (str,))

    structs_dir_list = match_struct_dirs(base_dir)

    if structs_names:
        structs_dir_list = list(filter(
            lambda path: os.path.basename(path) in structs_names,
            structs_dir_list
        ))

    nbr_structs = len(structs_dir_list)
    chunksize   = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)

    structs_data_list = list(
        filter(
            None,
            process_map(
                set_vasp_extractor,
                structs_dir_list,
                max_workers=workers,
                chunksize=chunksize,
                desc="Extracting infos from previous VASP output",
            ),
        )
    )

    structs_data = dict(structs_data_list)

    return structs_data


########################################


def vasp_output_structure(struct_dir: str) -> Tuple[Structure,Structure]|Tuple[None,None]:
    """
    Get a pymatgen Structure from the input and output of a VASP calculation.
    If the calculation did not converge, either by reaching max ionic steps, 
    by terminating on an error or by timeout, Tuple[None,None] is returned instead.

    Parameters:
        struct_dir (str):  Path to the calculation.

    Returns: (Structure, Structure)
        The initial and final structures if the calculation converged.
    """
    vasprun = _check_vasp_data(struct_dir)

    if vasprun is None:
        return (None, None)

    in_struct = vasprun.initial_structure
    out_struct = vasprun.final_structure

    return in_struct, out_struct


########################################


def batch_extract_vasp_structures(
    calc_dirs: List[str], workers: int = 1
) -> List[Tuple[Structure, Structure]]:
    """
    Get a list of structures from a list of path to VASP calculations.

    Parameters:
        calc_dirs (List[str]):  List of path to calculations.
        workers (int):  number of workers.

    Returns: (List[Tuple[Structure, Structure]])
        List of loaded structures.
    """
    check_type(workers, "workers", (int,))
    check_num_value(workers, "workers", ">", 0)

    return list(filter(
        lambda tup: tup != (None, None),
        process_map(
            vasp_output_structure,
            calc_dirs,
            max_workers=workers,
            desc="Extracting structures from VASP output",
        )
    ))


########################################


def dsol_calc_init(
        structure: Structure,
        calc_index: int,
        preset: Union[PMGStaticSet, Literal["DsolStaticSet"]] = "DSolStaticSet",
        user_corrections: Optional[Dict[str, Any]] = None,
    ) -> VaspInput:
    """
    Initializes one of the static calculations used for Δ-Sol method for one structure.

    Reference of the Δ-Sol method:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (values in Table I)

    Parameters:
        structure (Structure):      The input structure.

        calc_index (int):           An integer corresponding to a delta-sol static calculation:
                                    0 = E(N0), 
                                    1-2 = E(N0 + n), E(N0 - n) respectively, using N*_best, 
                                    3-4 = E(N0 + n), E(N0 - n) respectively, using N*_min, 
                                    5-6 = E(N0 + n), E(N0 - n) respectively, using N*_max.

        preset (str):               A pymatgen VASP static preset, or the homemade
                                    DSolStaticSet. Defaults to DSolStaticSet.

        user_corrections (dict):    Additional corrections provided by the user in a
                                    separate .yaml file.

    Returns:
        The corresponding VaspInput object.
    """
    check_type(structure, "structure", (Structure,))
    check_type(calc_index, "calc_index", (int,))
    check_num_value(calc_index, "calc_index", ">=", 0)
    check_num_value(calc_index, "calc_index", "<=", 6)
    if preset not in PMGStaticSet and preset != "DSolStaticSet":
        raise ValueError(
            "'preset' argument value is not a supported preset. "
            "Supported presets are:\n"
            f"{PMGStaticSet + set(('DSolStaticSet',))}\n"
            f"'preset' got value '{preset}' instead."
        )
    if user_corrections is not None:
        check_type(user_corrections, "user_corrections", (Dict,))

    nb_val_elec = get_all_valence_electrons(structure)
    run_set = vasp_static_settings(structure, preset, user_corrections)

    # Search for the right N* parameter to use with respect to the functional
    pot_func = run_set.get("POTCAR_FUNCTIONAL", "PBE")

    n_ratio = get_dsol_n_ratio(
        structure=structure,
        dft_functional=pot_func,
        n_star_idx=calc_index
    )

    nelect = nb_val_elec + n_ratio if calc_index % 2 == 1 else nb_val_elec - n_ratio

    if preset == "DSolStaticSet":
        run_set = vasp_static_settings(
            structure, preset, nelect=nelect, user_corrections=user_corrections
        )

    else:
        run_dict = run_set.as_dict()
        run_dict["INCAR"].update({"NELECT": nelect})
        run_set = VaspInput.from_dict(run_dict)

    return run_set


########################################
