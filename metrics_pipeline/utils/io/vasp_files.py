"""Management of VASP files for quantum simulations."""

import os
import re
from pathlib import Path
import typing as tp
from enum import Enum
import json
import warnings
import functools as ft

from monty.os import cd
import xml.etree.ElementTree as ET
from tqdm.contrib.concurrent import process_map

from pymatgen.core import SiteCollection, Structure
from pymatgen.io.vasp.sets import (
    VaspInput, VaspInputSet,
    MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPMetalRelaxSet, MVLRelax52Set, MVLScanRelaxSet,
    MPStaticSet, MPSOCSet, MatPESStaticSet, MPScanStaticSet
)
from pymatgen.io.vasp import VaspInput, Vasprun, Poscar, Xdatcar

from .io_base import PathLike, check_file_or_dir
# TODO: Split responsibilities better, e.g. define compound convenience functions
# outside the package.
from metrics_pipeline.utils.delta_sol import DSolStaticSet, get_dsol_n_ratio, get_all_valence_electrons


U_VALUES = {
    "F": {
        "Ag": 1.5, "Co": 3.4, "Cr": 3.5, "Cu": 4.0, #"Cu": 4 -> 4.0
        "Fe": 4.0, "Mn": 3.9, "Mo": 3.5, "Nb": 1.5, #"Mo": 4.38 -> 3.5 (according to the reference)
        "Ni": 6.0, "Re": 2.0, "Ta": 2.0, "V": 3.1,  #"Ni": 6 -> 6.0, "Re": 2 -> 2.0, "Ta": 2 -> 2.0
        "W": 4.0
    },
    "O": {
        "Ag": 1.5, "Co": 3.4, "Cr": 3.5, "Cu": 4.0, #"Cu": 4 -> 4.0
        "Fe": 4.0, "Mn": 3.9, "Mo": 3.5, "Nb": 1.5, #"Mo": 4.38 -> 3.5 (according to the reference)
        "Ni": 6.0, "Re": 2.0, "Ta": 2.0, "V": 3.1,  #"Ni": 6 -> 6.0, "Re": 2 -> 2.0, "Ta": 2 -> 2.0
        "W": 4.0                          
    },
    "S": {
        "Fe": 1.9, "Mn": 2.5
    }}

"""
Values of the Hubbard U correction used in GGA + U framework, as fitted by Jain et al.

Reference:
    A. Jain, G. Hautier, C.J. Moore, S.P. Ong, 
    C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder, 
    Computational Materials Science, 50, 2295-2310 (2011)
"""

class PresetEnum(Enum):
    """Base Enum with convenient methods to navigate between stored presets."""
    @classmethod
    def names(cls) -> list[str]:
        """List of all presets names."""
        return list(member.name for member in cls)

    @classmethod
    def values(cls) -> list[type[VaspInputSet]]:
        """List of all preset classes."""
        return list(member.value for member in cls)

    @classmethod
    def items(cls) -> list[tuple[str, type[VaspInputSet]]]:
        """List of name / preset pairs as tuples, similar to the dict.items() method."""
        return [(member.name, member.value) for member in cls]

    @classmethod
    def is_preset(cls, name: str) -> bool:
        """Whether given name corresponds to an existing preset (case insensitive)."""
        return name.upper() in [preset_name.upper() for preset_name in cls.names()]

    @classmethod
    def get_preset(cls, name: str) -> type[VaspInputSet]:
        """Get one of the stored callable preset classes by its name (case insensitive)."""
        if not cls.is_preset(name):
            raise ValueError(
                f"{name!r} is not a valid preset of {cls.__name__}. "
                f"Available presets in this class: {', '.join(cls.names())}."
            )
        for preset_name, preset in cls.items():
            if name.upper() == preset_name.upper():
                return preset

        raise RuntimeError(
            "Preset name was recognized by 'is_preset' but not found when looping. "
            "This error should never occur and is likely a bug. Please open an issue "
            "on the development repo of the library about this."
        )


class PMGRelaxSet(PresetEnum):
    """Enum class of known to date VASP relaxation presets implemented in pymatgen."""
    MITRELAXSET = MITRelaxSet
    MPRELAXSET = MPRelaxSet
    MPSCANRELAXSET = MPScanRelaxSet
    MPMETALRELAXSET = MPMetalRelaxSet
    MVLRELAX52SET = MVLRelax52Set
    MVLSCANRELAXSET = MVLScanRelaxSet


class PMGStaticSet(PresetEnum):
    """Enum class of known to date VASP static presets implemented in pymatgen."""
    MPSTATICSET = MPStaticSet
    MATPESSTATICSET = MatPESStaticSet
    MPSCANSTATICSET = MPScanStaticSet
    MPSOCSET = MPSOCSet


def _mitrelaxset_incar_corrections(n_sites: int|None = None) -> dict[str, tp.Any]:
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


def _relax_set_init(
    structure: SiteCollection,
    preset: PMGRelaxSet = PMGRelaxSet.MPRELAXSET,
    corrections: dict[str, tp.Any] | None = None,
) -> VaspInputSet:
    """Init a relaxation set of VASP input files."""
    corrections = {} if corrections is None else corrections
    incar_corrections = {}

    if preset == PMGRelaxSet.MITRELAXSET:
        incar_corrections = _mitrelaxset_incar_corrections(structure.num_sites)

    incar_corrections.update(corrections.get("INCAR", {}))

    vasp_input_set = PMGRelaxSet.get_preset(preset.name)(
        structure=structure,
        user_incar_settings = incar_corrections,
        user_kpoints_settings = corrections.get("KPOINTS", {}),
        user_potcar_settings = corrections.get("POTCAR", {}),
        user_potcar_functional = corrections.get("POTCAR_FUNCTIONAL", {}),
    )

    return vasp_input_set


def _static_set_init(
    structure: SiteCollection,
    preset: PMGStaticSet | tp.Literal["DSolStaticSet"] = PMGStaticSet.MPSTATICSET,
    nelect: float | None = None,
    corrections: dict[str, tp.Any] | None = None,
) -> VaspInputSet:
    """Init a static calculation set of VASP input files."""
    corrections = {} if corrections is None else corrections

    if preset == "DSolStaticSet":
        vasp_input_set = DSolStaticSet(
            structure=structure,
            incar_nelect=nelect,
            user_incar_settings=corrections.get("INCAR", {}),
            user_kpoints_settings=corrections.get("KPOINTS", {}),
            user_potcar_settings=corrections.get("POTCAR", {}),
            user_potcar_functional=corrections.get("POTCAR_FUNCTIONAL", {}),
        )
    else:
        vasp_input_set = PMGStaticSet.get_preset(preset.name)(
            structure=structure,
            user_incar_settings = corrections.get("INCAR", {}),
            user_kpoints_settings = corrections.get("KPOINTS", {}),
            user_potcar_settings = corrections.get("POTCAR", {}),
            user_potcar_functional = corrections.get("POTCAR_FUNCTIONAL", {}),
        )

    return vasp_input_set


def vasp_relaxation_settings(
    structure: SiteCollection,
    preset: PMGRelaxSet = PMGRelaxSet.MPRELAXSET,
    user_corrections: dict[str, tp.Any] | None = None,
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
    assert isinstance(structure, SiteCollection), TypeError(
        f"'structure' expected a type 'SiteCollection', got {type(structure).__name__}."
    )
    assert preset in PMGRelaxSet, (
        "'preset' argument not recognized. "
        "It must be one of the allowed pymatgen relaxation presets:\n"
        f"{', '.join(PMGRelaxSet.names())}."
    )

    vasp_input = _relax_set_init(
        structure, preset, user_corrections
    ).get_input_set()

    return vasp_input


def vasp_static_settings(
    structure: SiteCollection | None = None,
    preset: PMGStaticSet | tp.Literal["DSolStaticSet"] = PMGStaticSet.MPSTATICSET,
    nelect: float  |None = None,
    user_corrections: dict[str, tp.Any] | None = None,
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
        assert isinstance(structure, SiteCollection), TypeError(
            f"'structure' expected a type 'SiteCollection', got {type(structure).__name__}."
        )
    assert preset in PMGStaticSet or preset == "DSolStaticSet", (
        "'preset' argument not recognized. "
        "It must be 'DSolStaticSet' or one of the allowed pymatgen static presets:\n"
        f"{', '.join(PMGStaticSet.names())}."
    )

    vasp_input = _static_set_init(
        structure, preset, nelect, user_corrections
    ).get_input_set()

    return vasp_input


def _check_vasp_input(vasp_input: VaspInput) -> None:
    """Perform several tests to verify VaspInput correctness."""
    assert vasp_input.get("INCAR") is not None, ValueError(
            "check_vasp_input: There is no INCAR defined in the input !"
        )
    assert vasp_input.get("POSCAR") is not None, ValueError(
            "check_vasp_input: There is no POSCAR defined in the input !"
        )
    assert (
        vasp_input.get("KPOINTS") is not None or
        vasp_input["INCAR"].get("KSPACING") is not None
    ), ValueError("check_vasp_input: There is no KPOINTS or KSPACING tag defined in the input !")
    assert vasp_input.get("POTCAR") is not None, ValueError(
            "check_vasp_input: There is no POTCAR defined in the input !"
        )


def write_and_run_vasp(
    vasp_input: VaspInput, run_path: PathLike, vasp_exe: PathLike = "vasp"
) -> None:
    """
    Function to run VASP from a VaspInput object.

    Parameters
    ----------
    vasp_input: VaspInput
        The VaspInput object containing all necessary data to write VASP input files.

    run_path: str | Path
        Path to the directory where VASP files will be written and run.

    vasp_exe: str | Path
        Absolute path to the VASP executable. If not given, attempt the usual shortcut
        "vasp" launch command at given run_path.
    """
    _check_vasp_input(vasp_input)
    check_file_or_dir(run_path, "dir")
    if vasp_exe != "vasp":
        check_file_or_dir(vasp_exe, "file")

    vasp_input.write_input(output_dir=run_path)

    with cd(run_path):
        os.system(f"{vasp_exe}")


def get_struct_from_vasp(
    path: PathLike,
    try_vasprun: bool = True,
    try_contcar: bool = True,
    try_xdatcar: bool = True,
    try_poscar: bool = True
) -> Structure:
    """
    Attempt to extract the best structure possible from a VASP run directory.

    Algorithm
    ---------
    Attempt to extract a structure from a VASP run in following order,
    passing to the next step if the current one fails:
    - 1 - Tries to get the final structure from vasprun.xml
    - 2 - Tries to get the final structure from CONTCAR
    - 3 - Tries to get the last structure from XDATCAR
    - 4 - Tries to get the initial structure from POSCAR

    Parameters
    ----------
    path: str | Path
        Path to the VASP run directory.

    try_vasprun: bool
        Whether to try to get a structure from vasprun.xml. Defaults to True.

    try_contcar: bool
        Whether to try to get a structure from CONTCAR. Defaults to True.

    try_xdatcar: bool
        Whether to try to get a structure from XDATCAR. Defaults to True.

    try_poscar: bool
        Whether to try to get a structure from POSCAR. Defaults to True.

    Returns
    -------
    Structure
        A Structure object if it could be extracted.

    Raises
    ------
    `FileNotFoundError` if none of the attempts were successful.
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


def is_struct_dir(path: PathLike) -> bool:
    """
    Checks whether given path is a directory having naming convention of a structure directory.
    """
    #NOTE: Here 'path' must be a Path object, do not change for os.path.
    path = Path(path)
    return (
        path.is_dir()
        and re.match(r"\A[0-9]+_[A-Za-z0-9\(\)]+\Z", path.name) is not None
    )


def match_struct_dirs(
    path: PathLike,
    indices: list[int] | None = None,
    match_all: bool = True,
    unique: bool = True,
    no_return: bool = False
) -> list[str]:
    """
    Match structure directories at given path.
    
    Parameters
    ----------
    path: str | Path
        Where to search for structure directories.

    indices: list[int], optional
        Indices of wanted structure directories to match.
        If not given, matches all directories conforming to structure directory format.
        In this default behavior, match_all, unique and no_return are ignored.

    match_all: bool
        Whether each given index should match at least one directory.
        If True, raise an error when an index does not match any directory.
        Defaults to True.

    unique: bool
        Whether each given index should match at most one directory.
        If True, raise an error when an index matches several directories.
        Defaults to True.

    no_return: bool
        Whether to only match structure directories without returning the matching paths
        list for memory saving. If True, an empty list is always returned.
        Defaults to False.

    Returns
    -------
    list[str]
        List of paths of matching structure directories.

    Raises
    ------
    - `ValueError` if match_all is True and an index does not match any structure directory.
    - `ValueError` if unique is True and an index does match with several structure directories.
    """
    subtree = list(Path(path).iterdir())

    if indices is None:
        return sorted(list(map(str, filter(is_struct_dir, subtree))))

    all_matching_dirs = []

    for idx in sorted(indices):
        matching_dirs = list(
            filter(
                lambda f: f.name.startswith(f"{idx}_") and is_struct_dir(f),
                subtree
            )
        )

        if match_all and not matching_dirs:
            raise ValueError(
                f"Index '{idx}' do not match with any structure "
                f"directory in {path}."
            )
        if unique and len(matching_dirs) > 1:
            raise ValueError(
                f"Index '{idx}' matches with more than one structure directory, "
                "which should not be possible as structure indexation is unique "
                "inside a same batch. Please make sure that there are no parasite "
                f"directories in {path}, i.e. non-structure directories starting with "
                f"'{idx}_' or structure directories moved from other batches."
            )
        if not no_return:
            all_matching_dirs += sorted(matching_dirs)

    return list(map(str, all_matching_dirs))


def converged_vasprun(run_path: PathLike, **kwargs) -> Vasprun | None:
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


def _check_vasp_data(
    struct_dir: PathLike,
    path_to_summary: PathLike | None = None,
    key_to_check: str | None = None
) -> Vasprun | None:
    """
    Check whether given path leads to a previously accepted structure,
    its vasprun.xml file exists and is a normally terminated and converged run.

    Parameters
    ----------
    struct_dir: str | Path
        Directory containing a finished VASP calculation on a structure.

    path_to_summary: str
        Path to a JSON summary file containing results from previous step.
        Used to not consider structures that failed before.
        If not provided, all structure subdirs will be extracted.

    key_to_check: str
        The dict key associated to the bool used to verify eligibility.
        If path_to_summary is given, it must be given too.

    Returns
    -------
    Vasprun | None
        Vasprun object if it is eligible and normally converged, None otherwise.
    """
    check_file_or_dir(struct_dir, "dir")
    struct_dir = str(struct_dir)

    if path_to_summary is not None:
        if key_to_check is None:
            raise ValueError("If path_to_summary is given, key_to_check must be given too.")

        if not _check_summary_data(struct_dir, path_to_summary, key_to_check):
            return None

    return converged_vasprun(
        struct_dir, parse_dos=False, parse_eigen=False, parse_potcar_file=False
    )


def extract_vasp_data_for_convex_hull(
    struct_dir: PathLike,
    path_to_summary: PathLike | None = None,
    key_to_check: str | None = None
) -> tuple[str, dict[str, Structure | float]]:
    """
    Extracts VASP data from a previous run for one structure.
    Keeps only relevant data for relative stability calculation.

    Parameters
    ----------
    struct_dir: str | Path
        Directory containing a finished VASP calculation on a structure.
        
    path_to_summary: str
        Path to a JSON summary file containing results from previous step.
        Used to not consider structures that failed before.
        If not provided, all structure subdirs will be extracted.

    key_to_check: str
        The dict key associated to the bool used to verify eligibility.
        If path_to_summary is given, it must be given too.

    Returns
    -------
    tuple[str, dict]
        tuple containing the name of the struct_dir and corresponding dict,
        containing following data, used in stability calculation:
        - structure directory name,
        - structure chemical composition as Composition object,
        - generated structure energy (in eV).
    """
    vasprun = _check_vasp_data(struct_dir, path_to_summary, key_to_check)

    if vasprun is None:
        return "", {}

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


def extract_vasp_data_for_delta_sol_init(
    struct_dir: PathLike,
    path_to_summary: PathLike | None = None,
    key_to_check: str | None = None
) -> tuple[str, dict[str, Structure | float]]:
    """
    Extracts VASP data from a previous run for one structure.
    Keeps only relevant data for Δ-Sol method.

    Parameters
    ----------
    struct_dir: str | Path
        Directory containing a finished VASP calculation on a structure.
        
    path_to_summary: str
        Path to a JSON summary file containing results from previous step.
        Used to not consider structures that failed before.
        If not provided, all structure subdirs will be extracted.

    key_to_check: str
        The dict key associated to the bool used to verify eligibility.
        If path_to_summary is given, it must be given too.

    Returns
    -------
    tuple[str, dict]
        tuple containing the name of the struct_dir and corresponding dict,
        containing following data, used in stability calculation:
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


def extract_vasp_data_for_delta_sol_calc(
    struct_dir: PathLike,
    path_to_summary: PathLike | None = None,
    key_to_check: str | None = None
) -> tuple[str, dict[str, Structure | float]]:
    """
    Extract the results of Δ-Sol computations.
    
    Parameters
    ----------
    struct_dir: str | Path
        Directory containing a finished VASP calculation on a structure.
        
    path_to_summary: str
        Path to a JSON summary file containing results from previous step.
        Used to not consider structures that failed before.
        If not provided, all structure subdirs will be extracted.

    key_to_check: str
        The dict key associated to the bool used to verify eligibility.
        If path_to_summary is given, it must be given too.

    Returns
    -------
    tuple[str, dict]
        tuple containing the name of the struct_dir and corresponding dict,
        containing following data, used in stability calculation:
        - structure itself,
        - its final energy (in eV).
    """
    check_file_or_dir(struct_dir, "dir")

    calc_dirs   = list(filter(os.path.isdir, os.listdir(struct_dir)))
    struct_dict = {}

    for calc_dir in calc_dirs:
        if (
            path_to_summary is not None
            and not _check_summary_data(struct_dir, path_to_summary, key_to_check) # type: ignore
        ):
            return "", {}

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


def batch_extract_vasp_data(
    method: tp.Literal["convex_hull", "delta_sol_init", "delta_sol_calc"],
    base_dir: PathLike,
    structs_names: tp.Sequence[str] | None = None,
    path_to_summary: PathLike | None = None,
    key_to_check: str | None = None,
    workers: int | None = None
) -> dict[str, dict[str, tp.Any]]:
    """
    Extracts VASP data from a previous run for each structure directory in given directory.
    Keeps only data that are useful according to given method.

    Parameters
    ----------
    method: str
        Name of the method that will use the data, used to know which data should be extracted.
        urrently supported methods: 'convex_hull', 'delta_sol_init', 'delta_sol_calc'.

    base_dir: str | Path
        Directory containing structures subdirs to extract data from.

    structs_names: Sequence[str]
        Provide specific structures sub-directories to extract data from.
        If specified, only specified subdirs in base_dir are checked.
        If not, all subdirs in base_dir are checked.
        
    path_to_summary: str
        Path to a JSON summary file containing results from previous step.
        Used to not consider structures that failed before.
        If not provided, all structure subdirs will be extracted.

    key_to_check: str
        The dict key associated to the bool used to verify eligibility.
        If path_to_summary is given, it must be given too.

    workers: int, optional
        Number of parallel processes to spawn. If not given,
        `tqdm.contrib.concurrent.process_map()` default is used.
        Pass 0 to disable `process_map()` and execute sequentially.

    Returns
    -------
    dict[str, dict[str, Any]]
        Dict with structure directory names as keys, and a dict containing useful data
        according to chosen method for corresponding structure as values.
        See respective extraction functions for more details on extracted informations.
    """
    check_file_or_dir(base_dir, "dir")
    assert os.listdir(str(base_dir)), (
        f"{str(base_dir)}: Directory exists but is empty."
    )
    if workers is not None:
        assert isinstance(workers, int), TypeError(
            f"'workers' expected a type 'int', got {type(workers).__name__}."
        )
        assert workers >= 0, ValueError(
            f"'workers' must be positive or zero, got {workers}."
        )

    match method:
        case "convex_hull":
            extract_fn = extract_vasp_data_for_convex_hull
        case "delta_sol_init":
            extract_fn = extract_vasp_data_for_delta_sol_init
        case "delta_sol_calc":
            extract_fn = extract_vasp_data_for_delta_sol_calc
        case str():
            raise NotImplementedError(
                f"Provided method ({method}) is not supported.\n"
                "Supported methods are: "
                "'convex_hull', 'delta_sol_init', 'delta_sol_calc'."
            )
        case _:
            raise TypeError(f"'method' expected a type 'str', got {type(method).__name__}.")

    vasp_extractor = ft.partial(
        extract_fn, path_to_summary=path_to_summary, key_to_check=key_to_check
    )

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
                vasp_extractor,
                structs_dir_list,
                max_workers=workers,
                chunksize=chunksize,
                desc="Extracting infos from previous VASP output",
            ),
        )
    )
    structs_data = dict(structs_data_list)

    return structs_data


def vasp_output_structure(struct_dir: str) -> tuple[Structure,Structure]|tuple[None,None]:
    """
    Get a pymatgen Structure from the input and output of a VASP calculation.
    If the calculation did not converge, either by reaching max ionic steps, 
    by terminating on an error or by timeout, tuple[None,None] is returned instead.

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


def batch_extract_vasp_structures(
    calc_dirs: list[str], workers: int | None = None
) -> list[tuple[Structure, Structure]]:
    """
    Get a list of structures from a list of path to VASP calculations.

    Parameters:
        calc_dirs (list[str]):  List of path to calculations.

        workers (int):          Number of parallel processes to spawn.
                                If not given, tqdm.contrib.concurrent.process_map
                                default is used.

    Returns: (list[tuple[Structure, Structure]])
        list of loaded structures.
    """
    if workers is not None:
        assert isinstance(workers, int), TypeError(
            f"'workers' expected a type 'int', got {type(workers).__name__}."
        )
        assert workers >= 0, ValueError(
            f"'workers' must be positive or zero, got {workers}."
        )

    return list(filter(
        lambda tup: tup != (None, None),
        process_map(
            vasp_output_structure,
            calc_dirs,
            max_workers=workers,
            desc="Extracting structures from VASP output",
        )
    ))


def dsol_calc_init(
        structure: Structure,
        calc_index: int,
        preset: PMGStaticSet | tp.Literal["DSolStaticSet"] = "DSolStaticSet",
        user_corrections: dict[str, tp.Any] | None = None,
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
    assert isinstance(structure, Structure), TypeError(
        f"'structure' expected a type 'Structure', got {type(structure).__name__}."
    )
    assert isinstance(calc_index, int), TypeError(
        f"'calc_index' expected a type 'int', got {type(calc_index).__name__}."
    )
    assert 0 <= calc_index <= 6, ValueError(
        f"'calc_index' must be between 0 and 6 included."
    )
    if preset not in PMGStaticSet and preset != "DSolStaticSet":
        raise ValueError(
            f"'preset' got unsupported value {preset!r}. "
            f"Supported presets: {', '.join(PMGStaticSet.names() + ['DSolStaticSet'])}."
        )

    nb_val_elec = get_all_valence_electrons(structure)
    run_set = vasp_static_settings(structure, preset, user_corrections=user_corrections)

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