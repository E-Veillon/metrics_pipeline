"""
Functions to write VASP input files, launch VASP calculations and manage VASP output files.
"""

########################################
# SYSTEM I/O MODULES

import os
import json
import warnings
from copy import deepcopy
from typing import Optional, Dict, List, Union, Sequence, Tuple, Literal, Any
from pathlib import Path
from dataclasses import dataclass
import xml.etree.ElementTree as ET

########################################
# OPTIMIZATION MODULES

import re
#from scipy.constants import elementary_charge
from itertools import starmap, repeat#, chain
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure, SiteCollection
from pymatgen.core import SETTINGS
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import PotcarSingle
from pymatgen.io.vasp import VaspInput, Vasprun
from pymatgen.io.vasp.sets import (
    DictSet, MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPHSERelaxSet,
    MPMetalRelaxSet, MVLRelax52Set, MVLScanRelaxSet
)

########################################
# LOCAL MODULES

from screening_pipeline.utils.custom_types import (
    PathLike,
    PMGRelaxSetType,
    PMGStaticSetType,
    PMGRelaxSet,
    PMGStaticSet,
)
from screening_pipeline.utils.utils import _yaml_loader
from screening_pipeline.utils.matcher import flatten
from screening_pipeline.utils.fitted_values import U_VALUES
from screening_pipeline.utils.paths_io import batch_add_new_dirs
from screening_pipeline.utils.periodic_table import (
    get_all_valence_electrons,
    get_delta_sol_el_ratio,
)


########################################
# LOCAL CLASS


@dataclass
class DeltaSolStaticSet(MPRelaxSet):
    """
    Initialize VASP input files for Δ-Sol method computations using
    PBE_54_W_HASH pymatgen set of POTCAR files. Parameters are as 
    described in Δ-Sol method original work by Chan et al. in 2010.
    DFT+U corrections are used as proposed by Jain et al. in 2011.

    References:
        - M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
        (Ref 32 in screening_pipeline/Bibliography)

        - A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller, 
        K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).
        (Ref 14 in screening_pipeline/Bibliography)

    Args:
        structure (Structure):  The Structure to create inputs for. If None, the input
                                set is initialized without a Structure but one must be
                                set separately before the inputs are generated.

        incar_nelect (float):   The number of electrons to put in the NELECT INCAR tag.
                                In Δ-Sol, several computations with distinct number of 
                                electrons are done, this is a convenient arg to set that.
                                If not given, infers the Δ-Sol N0 electrons calculation 
                                from the given structure.

        **kwargs:               kwargs supported by DictSet.
    
    Raises: ValueError if neither structure nor nelect are given at instanciation time.
    """
    base_path = os.path.dirname(os.path.dirname(__file__))
    path = os.path.join(base_path, "config", "DeltaSolStaticSet.yaml")
    CONFIG = _yaml_loader(path, on_error="raise")

    def __init__(
            self,
            structure: Structure|None = None,
            incar_nelect: float|None = None,
            **kwargs
        ) -> None:
        """DeltaSolStaticSet init."""
        super().__init__(structure, **kwargs)

        if incar_nelect is None:
            try:
                incar_nelect = get_all_valence_electrons(structure)
            except TypeError as exc:
                raise ValueError("Either structure or incar_nelect must be given.") from exc

        self.incar_nelect = incar_nelect

    @property
    def incar_updates(self) -> Dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: Dict[str, Any] = {"MAGMOM": None, "NELECT": self.incar_nelect}
        return updates

if __name__ == "__main__":
    # Test for the DeltaSolStaticSet class. You may have to change the given path
    # to one pointing at a valid CIF structure file for it to work properly.
    # You'll also need to set PMG_VASP_PSP_DIR for POTCAR files in .pmgrc.yaml.
    test_path = "/home/elohan/screening-pipeline/screening_pipeline/_benchmarks/TiO2.cif"
    with open(test_path, "rt", encoding="utf-8") as test_file:
        struct = CifParser(test_file).parse_structures()[0]
    dset = DeltaSolStaticSet(struct).get_input_set()
    with open(
        os.path.join(DeltaSolStaticSet.base_path, "DeltaVaspInput.txt"),
        mode="wt", encoding="utf-8"
    ) as out:
        out.write(str(dset))
    print(dset.incar_nelect)
    dset.incar_nelect = 58
    print(dset.incar_nelect)

########################################
# LOCAL FUNCTIONS


def _get_potcar_enmax_values(
    structure: SiteCollection,
    default_config: Dict,
    user_potcar_dict: Optional[Dict] = None,
    user_potcar_functional: Optional[str] = None,
) -> List[float]:
    """
    Reminder:
    When installing the pipeline repository in a new environment, after having initialised
    POTCARs with pymatgen, assemble all ENMAX values in a file called 'ENMAX_<POTCAR_FUNCTIONAL>.txt'
    in the POTCARs pymatgen directory before calling this function.
    """

    potcar_dict = user_potcar_dict or default_config.get("POTCAR", {})
    functional = user_potcar_functional or default_config.get(
        "POTCAR_FUNCTIONAL", "PBE"
    )
    enmax_list = []

    enmax_path = os.path.join(
        SETTINGS["PMG_VASP_PSP_DIR"],
        PotcarSingle.functional_dir[functional],
        f"ENMAX_{functional}.txt",
    )

    if not Path(enmax_path).is_file():
        raise FileNotFoundError(
            f"{enmax_path} not found. Use gzip and grep commands "
            "to assemble ENMAX values in this file."
        )

    for elt in structure.composition.element_composition.elements:
        elt_symbol = str(elt)
        try:
            potcar_symbol = next(
                filter(lambda symbol: symbol.find(elt_symbol) != -1, potcar_dict.values())
            )
        except StopIteration as exc:
            raise ValueError(
                f"No POTCAR symbol found for element '{elt_symbol}' in provided POTCAR settings."
                "Please verify the .yaml file used for the calculation."
            ) from exc

        with open(enmax_path, "rt", encoding="utf-8") as enmax_file:

            try:
                enmax_line = next(
                    filter(
                        lambda line: f"POTCAR.{potcar_symbol}:" in line
                        or f"{potcar_symbol}/POTCAR:" in line,
                        enmax_file.readlines(),
                    )
                ).split()
            except StopIteration as exc:
                raise ValueError(
                    f"The POTCAR symbol '{potcar_symbol}' defined in the .yaml configuration file"
                    f"cannot be found in {enmax_path}."
                    "Make sure you are using the right POTCAR library or that the ENMAX file"
                    "compiles the right POTCARs."
                ) from exc

        enmax_value = enmax_line[enmax_line.index("ENMAX") + 2].rstrip(";")
        enmax_list.append(float(enmax_value))

    return enmax_list


########################################


def _MITRelaxSet_INCAR_corrections(**kwargs) -> Dict:
    """
    Corrects errors and imprecisions found in Pymatgen MITRelaxSet VASP preset's INCAR tags.
    Some tags should depend on the associated structure or POTCARs.

    parameters:
        structure (SiteCollection):     Structure associated with the VASP run.

        user_potcar_dict (dict):        If POTCAR corrections are overriding
                                        MITRelaxSet's defaults, they must be
                                        provided here to initialize right ENCUT tag value.
                                        If not provided, it will be initialized from original
                                        MITRelaxSet POTCARs.

        user_potcar_functional (str):   If POTCAR_FUNCTIONAL correction is overriding
                                        MITRelaxSet's default, it must be provided here
                                        to initialize right ENCUT tag value. If not provided,
                                        it will be initialized from MITRelaxSet default.

    Returns:
        A dictionnary containing the INCAR tags corrections for MITRelaxSet.
    """

    structure: SiteCollection = kwargs.pop("structure")
    default_config = MITRelaxSet.CONFIG
    user_potcar_dict: Dict = kwargs.pop("user_potcar_dict", None)
    user_potcar_functional: str = kwargs.pop("user_potcar_functional", None)
    enmax_list = _get_potcar_enmax_values(
        structure, default_config, user_potcar_dict, user_potcar_functional
    )

    corrected_ediff = round(float(5e-5) * structure.num_sites, 6)
    corrected_encut = 1.3 * max(enmax_list)
    corrected_ldaul = {
        "F": {
            "Ag": 2,
            "Co": 2,
            "Cr": 2,
            "Cu": 2,
            "Fe": 2,
            "Mn": 2,
            "Mo": 2,
            "Nb": 2,
            "Ni": 2,
            "Re": 2,
            "Ta": 2,
            "V": 2,
            "W": 2,
        },
        "O": {
            "Ag": 2,
            "Co": 2,
            "Cr": 2,
            "Cu": 2,
            "Fe": 2,
            "Mn": 2,
            "Mo": 2,
            "Nb": 2,
            "Ni": 2,
            "Re": 2,
            "Ta": 2,
            "V": 2,
            "W": 2,
        },
        "S": {
            "Fe": 2,
            "Mn": 2,  #"Mn": 2.5 -> 2 (quantum number l have to be an integer)
        },
    }
    corrected_ldauu = U_VALUES
    corrected_incar = {
        "EDIFF": corrected_ediff,
        "ENCUT": corrected_encut,
        "LDAUL": corrected_ldaul,
        "LDAUU": corrected_ldauu,
        "LMAXMIX": 4,  # Necessary to get reliable results with GGA + U framework on d-type orbitals
    }
    return corrected_incar


########################################

def vasp_launcher(vasp_exe: PathLike, path: PathLike, vasp_input: VaspInput) -> None:
    """
    Function to run VASP from a VaspInput object.

    Parameters:
        vasp_exe (str|Path):    Absolute path to the VASP executable.

        path (str|Path):        Path to the directory where VASP files will be written and run.

        vasp_input (VaspInput): The VaspInput object containing all necessary data to run VASP.
    """
    assert isinstance(vasp_exe, (Path, str))
    assert isinstance(path, (Path, str))
    assert isinstance(vasp_input, VaspInput)
    assert vasp_input.get("INCAR") is not None, (
    "vasp_launcher: There is no INCAR defined in the input !"
    )
    assert vasp_input.get("POSCAR") is not None, (
    "vasp_launcher: There is no POSCAR defined in the input !"
    )
    assert (
        vasp_input.get("KPOINTS") is not None
        or vasp_input["INCAR"].get("KSPACING") is not None
    ), "vasp_launcher: There is no KPOINTS or KSPACING tag defined in the input !"
    assert vasp_input.get("POTCAR") is not None, (
    "vasp_launcher: There is no POTCAR defined in the input !"
    )


    vasp_exe_list = list((vasp_exe,)) # Necessary for subprocess to take it as a full command
    path          = str(path)
    calc_dir      = path
    out_file      = os.path.join(path, "vasp.out")
    err_file      = os.path.join(path, "vasp.err")

    try:
        vasp_input.run_vasp(
            run_dir=calc_dir,
            vasp_cmd=vasp_exe_list,
            output_file=out_file,
            err_file=err_file,
        )
    except FileExistsError: # Case of several processors trying to do the same calculation at once
        print(f"WARNING: {str(path)}: This directory already exists, skipping...")

########################################

def _vasp_launcher_batch_wrapper(args_tuple: Tuple[PathLike, VaspInput]):
    vasp_exe   = args_tuple[0]
    path       = args_tuple[1]
    vasp_input = args_tuple[2]
    vasp_launcher(vasp_exe, path, vasp_input)


########################################


def vasp_batch_launch(
        vasp_exe: PathLike,
        base_dir: PathLike,
        inputs_data: Dict[PathLike, VaspInput],
        workers: int = 1
    ) -> None:
    """
    Creates a subdirectory with given names for each provided VaspInput object, 
    writes VASP input files in those subdirectories, then runs VASP inside each one.
    As this process can be very expensive as it actually does the VASP computations,
    it is recommended to parallelize it by setting workers > 1.

    Parameters:
        vasp_exe (Path|str):        Path to the VASP executable.

        inputs_data (dict):         Dict of structure data to run, of the form
                                    {subdir_name: VaspInput}.

        base_dir (str|Path):        The directory where the subdirs should be created.

        workers (int):              The number of parallel processes to spawn.
                                    Defaults to 1.
    """

    assert isinstance(
        vasp_exe, (Path, str)
    ), f"vasp_exe: expected a str or Path, got {type(vasp_exe)} instead."

    assert all(isinstance(vasp_input, VaspInput) for vasp_input in inputs_data.values()), \
    f"vasp_inputs: Expected VaspInput objects, got types listed below:\n \
    {print(list((type(vasp_input) for vasp_input in inputs_data.values())))}"

    base_dir = Path(base_dir)

    assert base_dir.is_dir()
    assert all(isinstance(subdir_name, (Path, str)) for subdir_name in inputs_data.keys())
    #assert len(vasp_inputs) == len(subdir_names)
    assert isinstance(workers, int) and workers >= 1

    subpaths_list = batch_add_new_dirs(base_dir=base_dir, new_subdirs=list(inputs_data.keys()))
    inputs_data   = {subpath: inputs_data.get(subpath.name) for subpath in subpaths_list}
    inputs_list   = list(zip(repeat(vasp_exe), inputs_data.items()))
    inputs_list   = list(tuple(flatten(input)) for input in inputs_list)

    process_map(
        _vasp_launcher_batch_wrapper,
        inputs_list,
        max_workers=workers,
        chunksize=1,
        desc="VASP computations",
    )


########################################


def _relax_set_init(
    structure: SiteCollection,
    preset: str = "MPRelaxSet",
    corrections: Optional[Dict] = None,
) -> DictSet:

    corrections = corrections or {}
    incar_corrections = {}

    if preset == "MITRelaxSet":
        incar_corrections = _MITRelaxSet_INCAR_corrections(
            structure=structure,
            user_potcar_dict=corrections.get("POTCAR") or None,
            user_potcar_functional=corrections.get("POTCAR_FUNCTIONAL") or None,
        )

    incar_corrections.update(corrections.get("INCAR", {}))
    kpoints_corrections = corrections.get("KPOINTS", {})
    potcar_corrections = corrections.get("POTCAR", {})
    potcar_functional_correction = corrections.get("POTCAR_FUNCTIONAL", {})

    if preset == "MITRelaxSet":
        return MITRelaxSet(
            structure=structure,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MPRelaxSet":
        return MPRelaxSet(
            structure=structure,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MPScanRelaxSet":
        return MPScanRelaxSet(
            structure=structure,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MPHSERelaxSet":
        return MPHSERelaxSet(
            structure=structure,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MPMetalRelaxSet":
        return MPMetalRelaxSet(
            structure=structure,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MVLRelax52Set":
        return MVLRelax52Set(
            structure=structure,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MVLScanRelaxSet":
        return MVLScanRelaxSet(
            structure=structure,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if isinstance(preset, str):
        raise ValueError(f"Provided string is not a valid preset name ({preset}).")
    
    else:
        raise TypeError(
            f"'preset' arg expected a str type, got {type(preset)} instead."
        )


########################################


def _static_set_init(
    struct_or_path: Union[SiteCollection, PathLike],
    from_prev_calc: bool = False,
    preset: str = "MPStaticSet",
    nelect: float|None = None,
    corrections: Optional[Dict] = None,
) -> DictSet:
    """"""
    from pymatgen.io.vasp.sets import MPStaticSet, MatPESStaticSet, MPScanStaticSet

    corrections = corrections or {}
    incar_corrections = corrections.get("INCAR", {})
    kpoints_corrections = corrections.get("KPOINTS", {})
    potcar_corrections = corrections.get("POTCAR", {})
    potcar_functional_correction = corrections.get("POTCAR_FUNCTIONAL", {})

    if from_prev_calc:
        dir_path = Path(struct_or_path)
        assert dir_path.is_dir()
        if preset == "DeltaSolStaticSet":
            return DeltaSolStaticSet.from_prev_calc(
                prev_calc_dir=dir_path,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        if preset == "MPStaticSet":
            return MPStaticSet.from_prev_calc(
                prev_calc_dir=dir_path,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        if preset == "MatPESStaticSet":
            return MatPESStaticSet.from_prev_calc(
                prev_calc_dir=dir_path,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        if preset == "MPScanStaticSet":
            return MPScanStaticSet.from_prev_calc(
                prev_calc_dir=dir_path,
                user_incar_settings=incar_corrections,
                user_kpoints_settings=kpoints_corrections,
                user_potcar_settings=potcar_corrections,
                user_potcar_functional=potcar_functional_correction,
            )
        elif isinstance(preset, str):
            raise ValueError(
                f"Provided string is not a valid preset name ({preset})."
            )
        else:
            raise TypeError(
                f"'preset' arg expected a str type, got {type(preset)} instead."
            )
    if preset == "DeltaSolStaticSet":
        return DeltaSolStaticSet(
            structure=struct_or_path,
            incar_nelect=nelect,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MPStaticSet":
        return MPStaticSet(
            structure=struct_or_path,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MatPESStaticSet":
        return MatPESStaticSet(
            structure=struct_or_path,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    if preset == "MPScanStaticSet":
        return MPScanStaticSet(
            structure=struct_or_path,
            user_incar_settings=incar_corrections,
            user_kpoints_settings=kpoints_corrections,
            user_potcar_settings=potcar_corrections,
            user_potcar_functional=potcar_functional_correction,
        )
    elif isinstance(preset, str):
        raise ValueError(f"Provided string is not a valid preset name ({preset}).")
    else:
        raise TypeError(
            f"'preset' arg expected a str type, got {type(preset)} instead."
        )


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
        "'structure' argument format not supported. "
        "It must be an instance of the SiteCollection class or one of its subclasses."
    )

    assert preset in PMGRelaxSet, (
        "'preset' argument not recognized. "
        "It must be one of the allowed pymatgen relaxation presets:\n"
        f"{PMGRelaxSet}"
    )

    assert (
        isinstance(user_corrections, dict) or user_corrections is None
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
    preset: PMGStaticSetType|"DeltaSolStaticSet" = "MPStaticSet",
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
                                        Can also be the homemade "DeltaSolStaticSet" if Δ-Sol
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

        nelect (float):                 Only useful if DeltaSolStaticSet is used.
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
    Check if given structure directory is present in the summary file and if it
    was rejected in the previous step.
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

def converged_vasprun(run_path: PathLike, **kwargs) -> Vasprun|None:
    """
    Read and parse vasprun.xml file at given VASP run directory.
    Check whether the VASP run terminated and converged correctly.
    Returns a Vasprun object if it is the case, or None otherwise.
    
    Parameters:
        run_path (Path|str):    Path to the directory containing the VASP calculation.

        **kwargs:               Additional keyword arguments to pass to pymatgen Vasprun class.
    
    Returns:
        Vasprun object if no error, timeout nor max ionic step reached was issued,
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
            f"{run_path}: No vasprun.xml file found at this location."
            "Make sure this file is present in its run directory."
        )

    try:
        vasprun = Vasprun(
            filename=vasprun_path,
            ionic_step_skip=kwargs.pop("ionic_step_skip", None),
            ionic_step_offset=kwargs.pop("ionic_step_offset", int(0)),
            parse_dos=kwargs.pop("parse_dos", True),
            parse_eigen=kwargs.pop("parse_eigen", True),
            parse_projected_eigen=kwargs.pop("parse_projected_eigen", False),
            parse_potcar_file=kwargs.pop("parse_potcar_file", True),
            occu_tol=kwargs.pop("occu_tol", float(1e-8)),
            separate_spins=kwargs.pop("separate_spins", False),
            exception_on_bad_xml=kwargs.pop("exception_on_bad_xml", True)
        )
    except ET.ParseError:
        return None
    except UnicodeDecodeError:
        warn_msg = f"WARNING: vasprun.xml file at {run_path} contains "
        warn_msg += "unreadable characters for 'utf-8' codec.\n"
        warn_msg += "Associated data is therefore considered erroneous and "
        warn_msg += "is not parsed further."
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
        #NOTE: Do not change type hint here, os.listdir() is not appropriate.
        return os.path.isdir(path) and re.match(
            r"\A[0-9]+_[A-Za-z0-9\(\)]+\Z", os.path.basename(path)
        ) is not None

    if method == "convex_hull":
        set_vasp_extractor = partial(
            extract_vasp_data_for_convex_hull,
            path_to_summary=path_to_summary
        )
    elif method == "delta_sol_init":
        set_vasp_extractor = partial(
            extract_vasp_data_for_delta_sol_init,
            path_to_summary=path_to_summary
        )
    elif method == "delta_sol_calc":
        set_vasp_extractor = partial(
            extract_vasp_data_for_delta_sol_calc,
            path_to_summary=path_to_summary
        )
    else:
        raise NotImplementedError(
            f"Provided method ({method}) is not supported.\n"
            "Supported methods are: 'convex_hull', 'delta_sol_init', 'delta_sol_calc'."
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


def delta_sol_inputs_init(
    structs_data: dict,
    preset: PMGStaticSetType|"DeltaSolStaticSet" = "DeltaSolStaticSet",
    user_corrections: Optional[Dict] = None,
    with_uncertainties: bool = False,
) -> Tuple[List]:
    """
    Initialize Vasp static input sets for structures with N0 - n and N0 + n electrons per cell,
    used in band gap calculations with Δ-Sol method.

    Reference:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (reference 32 in screening_pipeline/Bibliography, values in Table I)

    Parameters:
        structs_data (dict):        Dict containing relevant VASP data from a previous relaxation,
                                    as provided by the batch_extract_vasp_data function.

        preset (PMGStaticSet):      One of the pymatgen static VASP preset labeled 'StaticSet', or
                                    the homemade DeltaSolStaticSet. Defaults to DeltaSolStaticSet.

        user_corrections (dict):    User defined settings. It allows to override some of 
                                    the preset INCAR, KPOINTS or POTCAR settings if necessary.
                                    Defaults to None.

        with_uncertainties (bool):  Whether to also compute uncertainty boundaries of Δ-Sol method.
                                    Defaults to False.

    Returns:
        Tuple[List[str], List[VaspInput]]:
        A list of the subdirectory names corresponding to Δ-Sol runs,
        and a list of the VaspInput objects for these runs, both in same order.
    """

    def _inputs_init(name: str, data: Dict[str, Any]) -> List[Tuple[str, VaspInput]]:

        # Initialize input sets
        structure = data["structure"]
        N_val     = get_all_valence_electrons(structure)
        run_set   = vasp_static_settings(
            structure, preset=preset, user_corrections=user_corrections
        )

        # Search for the right N* parameter to use
        pot_func = run_set.get("POTCAR_FUNCTIONAL", "PBE")

        if "LDA" in pot_func:
            delta_sol_functional = "LDA"
        elif "PBE" in pot_func:
            delta_sol_functional = "PBE"
        elif "AM05" in pot_func:
            delta_sol_functional = "AM05"
        else:
            raise NotImplementedError(
                "Provided POTCAR functional is not implemented for Δ-Sol method.\n"
                "Supported functionals are 'LDA', 'PBE' and 'AM05'."
            )

        n_star_types = ("BEST", "MIN", "MAX") if with_uncertainties else ("BEST",)

        if preset != "DeltaSolStaticSet":
            run_N0_dict = run_set.as_dict()
            run_N0_dict["INCAR"].update({"NELECT": N_val})
            run_N0 = VaspInput.from_dict(run_N0_dict)
        else:
            run_N0 = deepcopy(run_set)

        run_N0_path = "_".join((name , "neutral"))
        struct_runs_list = [(run_N0_path, run_N0)]


        # Compute relevant n = N_val / N*
        for n_star_type in n_star_types:

            n_ratio = get_delta_sol_el_ratio(
                structure=structure,
                dft_functional=delta_sol_functional,
                n_star_type=n_star_type
            )
            if preset != "DeltaSolStaticSet":
                run_plus, run_minus = run_set.as_dict(), run_set.as_dict()

                run_plus["INCAR"].update({"NELECT": N_val + n_ratio})
                run_minus["INCAR"].update({"NELECT": N_val - n_ratio})

                run_plus  = VaspInput.from_dict(run_plus)
                run_minus = VaspInput.from_dict(run_minus)
            else:
                run_plus = vasp_static_settings(
                    structure=structure,
                    preset=preset,
                    nelect=N_val + n_ratio,
                    user_corrections=user_corrections
                )
                run_minus = vasp_static_settings(
                    structure=structure,
                    preset=preset,
                    nelect=N_val - n_ratio,
                    user_corrections=user_corrections
                )

            run_plus_path    = "_".join((name , n_star_type.lower(), "plus"))
            run_minus_path   = "_".join((name , n_star_type.lower(), "minus"))

            struct_runs_list += [(run_plus_path, run_plus), (run_minus_path, run_minus)]

        return struct_runs_list

    inputs_data = dict(flatten(list(starmap(_inputs_init, list(structs_data.items())))))

    return inputs_data


########################################


def _match_calc_index(calc_index: int) -> str|None:
    match calc_index:
        case 0:
            return None
        case 1|2:
            return "BEST"
        case 3|4:
            return "MIN"
        case 5|6:
            return "MAX"
        case int():
            raise ValueError("calc_index must be between 0 and 6 included.")
        case _:
            raise TypeError(f"Expected 'int' type, got '{type(calc_index)}' type instead.")

if __name__ == "__main__":
    # Unit test for _match_calc_index().
    print("Wanted output:")
    print("None\nBEST BEST\nMIN MIN\nMAX MAX")
    print("Actual output:")
    print(_match_calc_index(0))
    print(_match_calc_index(1), _match_calc_index(2))
    print(_match_calc_index(3), _match_calc_index(4))
    print(_match_calc_index(5), _match_calc_index(6))

########################################


def delta_sol_calculation_init(
        structure: Structure,
        calc_index: int,
        preset: str = "DeltaSolStaticSet",
        user_corrections: Optional[Dict[str, Any]] = None,
    ) -> VaspInput:
    """
    Initializes one of the static calculations used for delta-sol method for one structure.
    Parameters:
        structure (Structure):      The input structure.

        calc_index (int):           An integer corresponding to a delta-sol static calculation:
                                    0 = E(N0), 
                                    1-2 = E(N0 + n), E(N0 - n) respectively, using N*_best, 
                                    3-4 = E(N0 + n), E(N0 - n) respectively, using N*_min, 
                                    5-6 = E(N0 + n), E(N0 - n) respectively, using N*_max.

        preset (str):               A pymatgen VASP static preset, or the homemade
                                    DeltaSolStaticSet. Defaults to DeltaSolStaticSet.

        user_corrections (dict):    Additional corrections provided by the user in a
                                    separate .yaml file.
    
    Returns:
        The corresponding VaspInput object.
    """

    assert isinstance(structure, Structure)
    assert isinstance(calc_index, int) and (0 <= calc_index <= 6)
    assert preset in PMGStaticSet or preset == "DeltaSolStaticSet"
    assert isinstance(user_corrections, Dict) or user_corrections is None

    N_val = get_all_valence_electrons(structure)
    run_set = vasp_static_settings(structure, preset, user_corrections=user_corrections)

    # Search for the right N* parameter to use with respect to the functional
    pot_func = run_set.get("POTCAR_FUNCTIONAL", "PBE")

    if "LDA" in pot_func:
        delta_sol_functional = "LDA"
    elif "PBE" in pot_func:
        delta_sol_functional = "PBE"
    elif "AM05" in pot_func:
        delta_sol_functional = "AM05"
    else:
        raise NotImplementedError(
            "Provided POTCAR functional is not implemented for Δ-Sol method.\n"
            "Recognized functionals are 'LDA', 'PBE', and 'AM05'."
        )

    n_star_type = _match_calc_index(calc_index)

    if n_star_type is not None:
        n_ratio = get_delta_sol_el_ratio(
            structure=structure,
            dft_functional=delta_sol_functional,
            n_star_type=n_star_type
        )
        nelect = N_val + n_ratio if calc_index % 2 == 1 else N_val - n_ratio

        if preset == "DeltaSolStaticSet":
            run_set = vasp_static_settings(
                structure, preset, nelect=nelect, user_corrections=user_corrections
            )
    else:
        nelect = N_val

    if preset != "DeltaSolStaticSet":
        run_dict = run_set.as_dict()
        run_dict["INCAR"].update({"NELECT": nelect})
        run_set = VaspInput.from_dict(run_dict)

    return run_set
