"""
Functions to manage, read, and write VASP files and launch VASP computations.
"""


import os
import re
import json
import warnings
from itertools import repeat
from functools import partial
from typing import Optional, Dict, List, Union, Sequence, Tuple, Literal, Any
from pathlib import Path
from monty.os import cd
import xml.etree.ElementTree as ET
from tqdm.contrib.concurrent import process_map

# PYTHON MATERIAL GENOMICS
from pymatgen.core import Structure, SiteCollection
from pymatgen.io.vasp import VaspInput, Vasprun
from pymatgen.io.vasp.sets import (
    DictSet, MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPHSERelaxSet,
    MPMetalRelaxSet, MVLRelax52Set, MVLScanRelaxSet,
    MPStaticSet, MatPESStaticSet, MPScanStaticSet
)

# LOCAL IMPORTS
from custom_types import (
    PathLike, PMGRelaxSetType, PMGStaticSetType,
    PMGRelaxSet, PMGStaticSet,
)
from delta_sol import DSolStaticSet
from fitted_values import U_VALUES


########################################


def _mitrelaxset_incar_corrections(n_sites: int|None = None) -> Dict:
    """
    Systematic correction for MITRelaxSet INCAR tags that do not match with 
    parameters given in the original work of Jain et al.

    Reference:
        A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller,
        K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).
        (reference 14 in screening_pipeline/Bibliography)

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

    if not isinstance(vasp_input, VaspInput):
        raise TypeError(
            "'vasp_input' argument expected a type 'pymatgen.io.vasp.VaspInput', "
            f"got '{type(vasp_input)}' instead."
        )
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

    if not isinstance(run_path, (Path, str)):
        raise TypeError(
            "'path' argument expected a type 'str' or 'pathlib.Path', "
            f"got '{type(run_path)}' instead."
        )
    if not isinstance(vasp_exe, (Path, str)):
        raise TypeError(
            "'vasp_exe' argument expected a type 'str' or 'pathlib.Path', "
            f"got '{type(vasp_exe)}' instead."
        )

    vasp_input.write_input(output_dir=run_path)

    with cd(run_path):
        os.system(f"{vasp_exe}")


########################################


def _relax_set_init(
    structure: SiteCollection,
    preset: str = "MPRelaxSet",
    corrections: Optional[Dict] = None,
) -> DictSet:
    """Init a relaxation set of VASP input files."""
    corrections = corrections or {}
    incar_corrections = {}

    if preset == "MITRelaxSet":
        incar_corrections = _mitrelaxset_incar_corrections(structure.num_sites)

    incar_corrections.update(corrections.get("INCAR", {}))
    kpoints_corrections = corrections.get("KPOINTS", {})
    potcar_corrections = corrections.get("POTCAR", {})
    potcar_functional_correction = corrections.get("POTCAR_FUNCTIONAL", {})

    match preset:
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
        case str():
            raise ValueError(f"Provided string is not a valid preset name ({preset}).")
        case _:
            raise TypeError(
                f"'preset' arg expected a type 'str', got {type(preset)} instead."
            )
    return preset_obj

########################################


def _static_set_init(
    struct_or_path: Union[SiteCollection, PathLike],
    from_prev_calc: bool = False,
    preset: str = "MPStaticSet",
    nelect: float|None = None,
    corrections: Optional[Dict] = None,
) -> DictSet:
    """Init a static calculation set of VASP input files."""

    corrections = corrections or {}
    incar_corrections = corrections.get("INCAR", {})
    kpoints_corrections = corrections.get("KPOINTS", {})
    potcar_corrections = corrections.get("POTCAR", {})
    potcar_functional_correction = corrections.get("POTCAR_FUNCTIONAL", {})

    if from_prev_calc:
        dir_path = str(struct_or_path)
        assert os.path.isdir(dir_path)

        match preset:
            case "DSolStaticSet":
                preset_obj = DSolStaticSet.from_prev_calc(
                    prev_calc_dir=dir_path,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )
            case "MPStaticSet":
                preset_obj = MPStaticSet.from_prev_calc(
                    prev_calc_dir=dir_path,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )
            case "MatPESStaticSet":
                preset_obj = MatPESStaticSet.from_prev_calc(
                    prev_calc_dir=dir_path,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )
            case "MPScanStaticSet":
                preset_obj = MPScanStaticSet.from_prev_calc(
                    prev_calc_dir=dir_path,
                    user_incar_settings=incar_corrections,
                    user_kpoints_settings=kpoints_corrections,
                    user_potcar_settings=potcar_corrections,
                    user_potcar_functional=potcar_functional_correction,
                )
            case str():
                raise ValueError(
                    f"Provided string is not a valid preset name ({preset})."
                )
            case _:
                raise TypeError(
                    f"'preset' arg expected a str type, got {type(preset)} instead."
                )
        return preset_obj

    match preset:
        case "DSolStaticSet":
            preset_obj = DSolStaticSet(
                structure=struct_or_path,
                incar_nelect=nelect,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MPStaticSet":
            preset_obj = MPStaticSet(
                structure=struct_or_path,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MatPESStaticSet":
            preset_obj = MatPESStaticSet(
                structure=struct_or_path,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case "MPScanStaticSet":
            preset_obj = MPScanStaticSet(
                structure=struct_or_path,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        case str():
            raise ValueError(f"Provided string is not a valid preset name ({preset}).")
        case _:
            raise TypeError(
                f"'preset' arg expected a str type, got {type(preset)} instead."
            )
    return preset_obj

########################################


def vasp_relaxation_settings(
    structure: SiteCollection,
    preset: PMGRelaxSetType = "MPRelaxSet",
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
    assert isinstance(structure, SiteCollection), (
        "'structure' argument expected a type "
        "'pymatgen.core.structure.SiteCollection', "
        f"got '{type(structure)}' instead."
    )
    assert preset in PMGRelaxSet, (
        "'preset' argument not recognized. "
        "It must be one of the allowed pymatgen relaxation presets:\n"
        f"{PMGRelaxSet}"
    )
    assert (
        isinstance(user_corrections, dict)
        or user_corrections is None
    ), "user_incar_settings must be a dict or None"

    vasp_input = _relax_set_init(
        structure=structure,
        preset=preset,
        corrections=user_corrections
    ).get_input_set()

    return vasp_input


########################################


def vasp_static_settings(
    structure: Optional[SiteCollection] = None,
    preset: PMGStaticSetType|"DSolStaticSet" = "MPStaticSet",
    from_prev_calc: bool = False,
    prev_calc_dir: Optional[PathLike] = None,
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

        from_prev_calc (bool):          Whether to get final structure, INCAR and KPOINTS settings
                                        from a previous VASP run. INCAR tags will still be managed
                                        to fit a static calculation if previous run is a relaxation.
                                        For the sake of consistency, it is recommended to use the
                                        static preset corresponding to previous relaxation preset
                                        in this case (e.g. MPStaticSet for a relaxation with
                                        MPRelaxSet). If set to True, a directory to extract data
                                        from must be provided and structure argument is ignored.
                                        Defaults to False.

        prev_calc_dir (str|Path):       Directory to extract previous VASP run data from when 
                                        from_prev_calc is True. Otherwise, this argument is
                                        ignored.

        nelect (float):                 Only useful if DSolStaticSet is used.
                                        Sets the NELECT tag in INCAR file.
                                        Ignored if from_prev_calc is True.

        user_corrections (dict):        User defined settings. It allows to override some of 
                                        the preset INCAR, KPOINTS or POTCAR settings if necessary.
                                        Defaults to None.
    """
    assert isinstance(structure, SiteCollection) or structure is None, (
        "'structure' argument format not supported. "
        "It must be an instance of the SiteCollection "
        "class or one of its subclasses."
    )

    assert preset in PMGStaticSet, (
        "'preset' argument not recognized. "
        "It must be one of the allowed pymatgen static presets:\n"
        f"{PMGStaticSet}."
    )

    assert isinstance(user_corrections, dict) or user_corrections is None, (
        "user_corrections must be a dict or None"
    )

    if not from_prev_calc:
        vasp_input = _static_set_init(
            struct_or_path=structure,
            preset=preset,
            nelect=nelect,
            corrections=user_corrections
        ).get_input_set()

    else:
        assert isinstance(prev_calc_dir, (Path, str)), (
            "'from_prev_calc' was set to True, prev_calc_dir must be provided "
            "as 'str' or 'Path' type."
        )

        prev_calc_dir = Path(prev_calc_dir)

        assert (
            prev_calc_dir.is_dir()
        ), f"Prev_calc_dir: {prev_calc_dir} is not a valid directory."

        vasp_input = _static_set_init(
            struct_or_path=prev_calc_dir,
            from_prev_calc=from_prev_calc,
            preset=preset,
            corrections=user_corrections
        ).get_input_set()

    return vasp_input


########################################


def _check_summary_data(struct_dir: PathLike, summary_file: PathLike) -> bool:
    """
    Check whether given structure directory is present in the summary file and
    whether it was rejected in the previous step.
    """

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

    is_good_struct = (
        "True" in prev_struct_data.values() or
        "true" in prev_struct_data.values()
    )
    return is_good_struct


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
    if not isinstance(run_path, (Path,str)):
        raise TypeError(
            "Given 'struct-dir' argument value expected 'Path' or 'str', "
            f"got '{type(run_path)}' instead."
        )
    if not os.path.isdir(run_path):
        raise ValueError(f"{run_path}: No such directory found.")

    vasprun_path = os.path.join(run_path, "vasprun.xml")

    if not os.path.isfile(vasprun_path):
        raise ValueError(
            f"{run_path}: No vasprun.xml file found at this location. "
            "Make sure this file is present in its run directory."
        )

    try:
        vasprun = Vasprun(filename=vasprun_path, **kwargs)
    except ET.ParseError:
        return None
    except UnicodeDecodeError:
        warn_msg = (
            f"WARNING: vasprun.xml file at {run_path} contains "
            "unreadable characters for 'utf-8' codec.\n"
            "Associated data is therefore considered erroneous "
            "and is not parsed further."
        )
        warnings.warn(warn_msg)
        return None

    if not vasprun.converged:
        return None

    return vasprun


########################################


def extract_vasp_data_for_convex_hull(
    struct_dir: PathLike = ".",
    path_to_summary: Optional[PathLike] = None
) -> Tuple[str, Dict[str, Union[Structure, float]]]:
    """
    Extracts VASP data from a previous run for one structure.
    Keeps only relevant data for relative stability calculation.

    Parameters:
        struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.
        
        path_to_summary (str):  Path to a JSON summary file containing results from previous step.
                                Used to not consider structures that failed before.
                                If not provided, all structure subdirs will be extracted.

    Returns:
        Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict,
                                containing following data, used in stability calculation:
                                    - structure directory name,
                                    - structure chemical composition as Composition object,
                                    - generated structure energy (in eV),
    """

    assert isinstance(struct_dir, (Path, str))
    assert os.path.isdir(str(struct_dir))

    struct_dir: str = str(struct_dir)

    if path_to_summary is not None and not _check_summary_data(struct_dir, path_to_summary):
        return {}

    vasprun = converged_vasprun(
        struct_dir, parse_dos=False, parse_eigen=False, parse_potcar_file=False
    )

    if vasprun is None:
        return {}

    # We want data of the generated structure for convex hulls, not the relaxed one
    struct_name  = os.path.basename(struct_dir)
    generated_energy: float = vasprun.ionic_steps[0]["e_0_energy"]
    composition = vasprun.initial_structure.composition

    struct_dict = {
        "entry_id": struct_name,
        "composition": composition,
        "final_energy": generated_energy
    }
    struct_data = (struct_name, struct_dict)

    return struct_data


########################################


def extract_vasp_data_for_delta_sol_init(
    struct_dir: PathLike = ".",
    path_to_summary: Optional[PathLike] = None
) -> Tuple[str, Dict[str, Union[Structure, float]]]:
    """
    Extracts VASP data from a previous run for one structure.
    Keeps only relevant data for Δ-Sol method.

    Parameters:
        struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.
        
        path_to_summary (str):  Path to a JSON summary file containing results from previous step.
                                Used to not consider structures that failed before.
                                If not provided, all structure subdirs will be extracted.


    Returns:
        Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict,
                                containing following data, used for Δ-Sol method:
                                    - structure itself,
                                    - its final energy (in eV).
    """

    assert isinstance(struct_dir, (Path, str)), (
        TypeError(
            "'struct_dir' argument expected 'Path' or 'str' type, "
            "got '{type(struct_dir)}' instead."
        )
    )
    assert os.path.isdir(str(struct_dir)), (
        ValueError(f"{struct_dir}: No such directory found.")
    )

    struct_dir: Path = Path(struct_dir)

    if path_to_summary is not None and not _check_summary_data(struct_dir, path_to_summary):
        return {}

    vasprun = converged_vasprun(
        struct_dir, parse_dos=False, parse_eigen=False, parse_potcar_file=False
    )
    if vasprun is None:
        return {}

    struct_name = struct_dir.name
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
    struct_dir: PathLike = ".",
    path_to_summary: Optional[PathLike] = None
) -> Tuple[str, Dict[str, Union[Structure, float]]]:
    """
    Extract the results of Δ-Sol computations.
    
    Parameters:
        struct_dir (str|Path):  Structure directory containing subdirs of Δ-Sol computations.
        
        path_to_summary (str):  Path to a JSON summary file containing results from previous step.
                                Used to not consider structures that failed before.
                                If not provided, all structure subdirs will be extracted.


    Returns:
        Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict,
                                containing following data, used for Δ-Sol method:
                                    - structure itself,
                                    - its final energy (in eV).
    """
    assert isinstance(struct_dir, (Path, str))
    assert os.path.isdir(str(struct_dir))

    calc_dirs   = list(filter(lambda path: os.path.isdir(path), os.listdir(struct_dir)))
    struct_dict = {}

    for calc_dir in calc_dirs:
        if path_to_summary is not None and not _check_summary_data(struct_dir, path_to_summary):
            return {}

        calc_data = extract_vasp_data_for_delta_sol_init(
            struct_dir=calc_dir, path_to_summary=path_to_summary
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
        base_dir: PathLike = ".",
        structs_names: Optional[Sequence[str]] = None,
        path_to_summary: Optional[PathLike] = None,
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

    assert isinstance(base_dir, (Path, str))
    assert os.path.isdir(str(base_dir))
    assert os.listdir(str(base_dir))
    assert isinstance(path_to_summary, (Path, str)) or path_to_summary is None
    assert isinstance(workers, int) and workers >= 1

    def is_struct_dir(path: Path) -> bool:
        #NOTE: Do not change type hint to str here, os.listdir() is not appropriate.
        return os.path.isdir(path) and re.match(
            r"\A[0-9]+_[A-Za-z0-9\(\)]+\Z", os.path.basename(path)
        ) is not None

    match method:
        case "convex_hull":
            set_vasp_extractor = partial(
                extract_vasp_data_for_convex_hull,
                path_to_summary=path_to_summary
            )
        case "delta_sol_init":
            set_vasp_extractor = partial(
                extract_vasp_data_for_delta_sol_init,
                path_to_summary=path_to_summary
            )
        case "delta_sol_calc":
            set_vasp_extractor = partial(
                extract_vasp_data_for_delta_sol_calc,
                path_to_summary=path_to_summary
            )
        case str():
            raise NotImplementedError(
                f"Provided method ({method}) is not supported.\n"
                "Supported methods are: "
                "'convex_hull', 'delta_sol_init', 'delta_sol_calc'."
            )
        case _:
            raise TypeError(
                f"'method' expected a type 'str', got '{type(method)}' instead."
            )

    #NOTE: Do not change the Path object here, os.listdir() is not appropriate
    structs_dir_list = list(filter(is_struct_dir, Path(base_dir).iterdir()))

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

    if not isinstance(struct_dir, str):
        raise TypeError(
            f"'struct_dir' arg expected a 'str', got '{type(struct_dir)}' instead."
        )
    if not os.path.isdir(struct_dir):
        raise ValueError(
            f"{struct_dir}: No such directory found."
        )

    vasprun = converged_vasprun(
        struct_dir, parse_dos=False, parse_eigen=False, parse_potcar_file=False
    )
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
