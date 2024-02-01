'''
Functions to write VASP input files, launch VASP calculations and manage VASP output files.
'''


########################################
# SYSTEM I/O MODULES

import os
from monty.os.path import zpath
from typing import Optional, Dict, List, Union, Sequence, Tuple, Literal, Any
from pathlib import Path

########################################
# OPTIMIZATION MODULES

from scipy.constants import elementary_charge
from itertools import chain, starmap, repeat
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure, SiteCollection, Composition
from pymatgen.core import SETTINGS
from pymatgen.io.vasp.inputs import PotcarSingle
from pymatgen.io.vasp import VaspInput
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Chgcar, Oszicar
from pymatgen.io.vasp.sets import DictSet, MITRelaxSet

########################################
# LOCAL MODULES

from screening_pipeline.utils.utils import is_float
from screening_pipeline.utils.typing import PathLike, PMGRelaxSetType, PMGStaticSetType, PMGRelaxSet, PMGStaticSet
from screening_pipeline.utils.matcher import flatten
from screening_pipeline.utils.fitted_values import U_VALUES
from screening_pipeline.utils.paths import batch_add_new_dirs
from screening_pipeline.utils.periodic_table import get_all_valence_electrons, get_delta_sol_el_ratio

########################################
# LOCAL FUNCTIONS

def _get_POTCAR_ENMAX_values(
        structure: SiteCollection, 
        default_config: Dict, 
        user_potcar_dict: Optional[Dict] = None, 
        user_potcar_functional: Optional[str] = None
) -> List[float]:
    '''
    Reminder:
    When installing the pipeline repository in a new environment, after having initialised 
    POTCARs with pymatgen, assemble all ENMAX values in a file called 'ENMAX_<POTCAR_FUNCTIONAL>.txt 
    in the POTCARs pymatgen directory before calling this function.
    '''

    potcar_dict    = user_potcar_dict or default_config.get('POTCAR', {})
    functional     = user_potcar_functional or default_config.get('POTCAR_FUNCTIONAL', 'PBE')
    ENMAX_list     = []

    enmax_path = os.path.join(
        SETTINGS['PMG_VASP_PSP_DIR'], 
        PotcarSingle.functional_dir[functional], 
        f'ENMAX_{functional}.txt'
    )

    if not Path(enmax_path).is_file():
        raise FileNotFoundError(
            f'{enmax_path} not found. Use gzip and grep commands to assemble ENMAX values in this file.'
        )

    for elt in structure.composition.element_composition.elements:
        
        try:
            potcar_symbol = next(filter(
                lambda symbol: symbol.find(str(elt)) != -1, 
                potcar_dict.values()
            ))
            print(f'potcar symbol found: {potcar_symbol}')
        except StopIteration:
            raise ValueError(
                f"No POTCAR symbol found for element '{elt}' in provided POTCAR settings."
                "Please verify the .yaml file used for the calculation."
            )

        with open(enmax_path, 'r') as enmax_file:

            try:
                enmax_line = next(filter(
                    lambda line: f'POTCAR.{potcar_symbol}:' in line or f'{potcar_symbol}/POTCAR:' in line, 
                    enmax_file.readlines()
                )).split()
                print(f'enmax_line found: {enmax_line}')
            except StopIteration:
                raise ValueError(
                    f"The POTCAR symbol '{potcar_symbol}' defined in the .yaml configuration file"
                    f"cannot be found in {enmax_path}."
                    "Make sure you are using the right POTCAR library or that the ENMAX file"
                    "compiles the right POTCARs."
                )

        enmax_value = enmax_line[enmax_line.index('ENMAX') + 2].rstrip(';')
        print(f"Final ENMAX value found for '{potcar_symbol}': {enmax_value}")
        ENMAX_list.append(float(enmax_value))

    return ENMAX_list

########################################

def _MITRelaxSet_INCAR_corrections(**kwargs) -> Dict:
    '''
    Corrects errors and imprecisions found in Pymatgen MITRelaxSet VASP preset's INCAR tags.
    Some tags should depend on the associated structure or POTCARs.

    parameters:
        structure (SiteCollection):     Structure associated with the VASP run.

        user_potcar_dict (dict):        If POTCAR corrections are overriding MITRelaxSet's defaults, 
                                        they must be provided here to initialize right ENCUT tag value.
                                        If not provided, it will be initialized from original MITRelaxSet
                                        POTCARs.

        user_potcar_functional (str):   If POTCAR_FUNCTIONAL correction is overriding MITRelaxSet's default, 
                                        it must be provided here to initialize right ENCUT tag value.
                                        If not provided, it will be initialized from MITRelaxSet default.

    Returns:
        A dictionnary containing the INCAR tags corrections for MITRelaxSet.
    '''

    structure: SiteCollection   = kwargs.pop('structure')
    default_config              = MITRelaxSet.CONFIG
    user_potcar_dict: Dict      = kwargs.pop('user_potcar_dict', None)
    user_potcar_functional: str = kwargs.pop('user_potcar_functional', None)
    ENMAX_list = _get_POTCAR_ENMAX_values(
        structure, default_config, user_potcar_dict, user_potcar_functional
    )

    corrected_EDIFF = float(5e-5)*structure.num_sites
    corrected_ENCUT = 1.3*max(ENMAX_list)
    print(f"Calculated ENCUT: {corrected_ENCUT}")
    corrected_LDAUL = {
        'F': {
            'Ag': 2, 'Co': 2, 'Cr': 2, 'Cu': 2, 'Fe': 2, 
            'Mn': 2, 'Mo': 2, 'Nb': 2, 'Ni': 2, 'Re': 2, 
            'Ta': 2, 'V': 2, 'W': 2
        }, 
        'O': {
            'Ag': 2, 'Co': 2, 'Cr': 2, 'Cu': 2, 'Fe': 2, 
            'Mn': 2, 'Mo': 2, 'Nb': 2, 'Ni': 2, 'Re': 2, 
            'Ta': 2, 'V': 2, 'W': 2
        }, 
        'S': {
            'Fe': 2, 'Mn': 2 #'Mn': 2.5 -> 2 (quantum number l have to be an integer)
        }}
    corrected_LDAUU = U_VALUES
    corrected_INCAR = {
        "EDIFF": corrected_EDIFF,
        "ENCUT": corrected_ENCUT,
        "LDAUL": corrected_LDAUL,
        "LDAUU": corrected_LDAUU, 
        "LMAXMIX": 4 #Necessary to get reliable results with GGA + U framework on d-type orbitals
        }
    return corrected_INCAR

########################################

def vasp_launcher(vasp_exe: PathLike, vasp_input: VaspInput, path: PathLike) -> None:

    assert isinstance(vasp_exe, (Path, str))
    assert isinstance(vasp_input, VaspInput)
    assert isinstance(path, PathLike)

    vasp_exe_list = list((vasp_exe,))
    path          = str(path)
    calc_dir      = Path(path)
    out_file      = Path('/'.join((path, "vasp.out")))
    err_file      = Path('/'.join((path, "vasp.err")))
    
    try:
        vasp_input.run_vasp(
            run_dir=calc_dir, 
            vasp_cmd=vasp_exe_list, 
            output_file=out_file, 
            err_file=err_file
        )
    except FileExistsError:
        print(f'WARNING: {str(path)}: This directory already exists, skipping...')

########################################

def _vasp_launcher_batch_wrapper(args_tuple: Tuple[VaspInput, PathLike]):
    vasp_exe   = args_tuple[0]
    vasp_input = args_tuple[1]
    path       = args_tuple[2]
    vasp_launcher(vasp_exe, vasp_input, path)

########################################

def vasp_batch_launch(
        vasp_exe: PathLike, 
        vasp_inputs: Sequence[VaspInput], 
        base_dir: PathLike, 
        subdir_names: Sequence[PathLike], 
        workers: int = 1
    ) -> None:
    '''
    Creates a subdirectory with given names for each provided VaspInput object, 
    writes VASP input files in those subdirectories, then runs VASP inside each one.
    As this process can be very expensive as it actually does the VASP computations, 
    it is recommended to parallelize it by setting workers > 1.

    Parameters:
        vasp_exe (Path|str):        Path to the VASP executable.

        vasp_inputs ([VaspInput]):  The objects defining how to write VASP input files in each subdirectory.

        base_dir (str|Path):        The directory where the subdirs should be created.

        subdir_names ([str|Path]):  The names or subpaths for created subdirectories.
                                    Note that its length must match the length of vasp_inputs.

        workers (int):              The number of parallel processes to spawn.
                                    Defaults to 1.
    '''

    assert isinstance(vasp_exe, (Path, str))
    assert all(isinstance(vasp_input, VaspInput) for vasp_input in vasp_inputs)
    assert isinstance(base_dir, PathLike)
    
    base_dir = Path(base_dir)

    assert base_dir.is_dir()
    assert all(isinstance(subdir_name, PathLike) for subdir_name in subdir_names)
    assert len(vasp_inputs) == len(subdir_names)
    assert isinstance(workers, int) and workers >= 1

    subpaths_list = batch_add_new_dirs(base_dir=base_dir, new_subdirs=subdir_names)
    inputs_list   = list(zip(repeat(vasp_exe), vasp_inputs, subpaths_list))

    process_map(
        _vasp_launcher_batch_wrapper, 
        inputs_list, 
        max_workers=workers, 
        chunksize=1, 
        desc='VASP computations'
    )

########################################

def _RelaxSet_init(
        structure: SiteCollection, 
        preset: str = 'MITRelaxSet', 
        corrections: Optional[Dict] = None
    ) -> DictSet:

    from pymatgen.io.vasp.sets import   MITRelaxSet, MPRelaxSet, MPScanRelaxSet, \
                                        MPHSERelaxSet, MPMetalRelaxSet, MVLRelax52Set, \
                                        MVLScanRelaxSet

    corrections = corrections or {}
    incar_corrections = {}

    if preset == 'MITRelaxSet':
        incar_corrections = _MITRelaxSet_INCAR_corrections(
            structure=structure, 
            user_potcar_dict=corrections.get('POTCAR') or None, 
            user_potcar_functional=corrections.get('POTCAR_FUNCTIONAL') or None
        )

    incar_corrections.update(corrections.get('INCAR', {}))
    kpoints_corrections = corrections.get('KPOINTS', {})
    potcar_corrections  = corrections.get('POTCAR', {})
    potcar_functional_correction = corrections.get('POTCAR_FUNCTIONAL', {})

    match preset:
        case 'MITRelaxSet': return MITRelaxSet(
            structure=structure, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MPRelaxSet': return MPRelaxSet(
            structure=structure, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MPScanRelaxSet': return MPScanRelaxSet(
            structure=structure, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MPHSERelaxSet': return MPHSERelaxSet(
            structure=structure, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MPMetalRelaxSet': return MPMetalRelaxSet(
            structure=structure, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MVLRelax52Set': return MVLRelax52Set(
            structure=structure, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MVLScanRelaxSet': return MVLScanRelaxSet(
            structure=structure, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case str(): raise ValueError(f'Provided string is not a valid preset name ({preset}).')
        case _: raise TypeError(f'"preset" arg expected a str type, got {type(preset)} instead.')

########################################

def _StaticSet_init(
        struct_or_path: Union[SiteCollection, PathLike], 
        from_prev_calc: bool = False, 
        preset: str = 'MPStaticSet', 
        corrections: Optional[Dict] = None
    ) -> DictSet:

    from pymatgen.io.vasp.sets import MPStaticSet, MatPESStaticSet, MPScanStaticSet

    corrections = corrections or {}
    incar_corrections   = corrections.get('INCAR', {})
    kpoints_corrections = corrections.get('KPOINTS', {})
    potcar_corrections  = corrections.get('POTCAR', {})
    potcar_functional_correction = corrections.get('POTCAR_FUNCTIONAL', {})

    if from_prev_calc:
        dir_path = Path(struct_or_path)
        assert dir_path.is_dir()

        match preset:
            case 'MPStaticSet': return MPStaticSet.from_prev_calc(
                prev_calc_dir=dir_path, 
                user_incar_settings=incar_corrections, 
                user_kpoints_settings=kpoints_corrections, 
                user_potcar_settings=potcar_corrections, 
                user_potcar_functional=potcar_functional_correction
            )
            case 'MatPESStaticSet': return MatPESStaticSet.from_prev_calc(
                prev_calc_dir=dir_path, 
                user_incar_settings=incar_corrections, 
                user_kpoints_settings=kpoints_corrections, 
                user_potcar_settings=potcar_corrections, 
                user_potcar_functional=potcar_functional_correction
            )
            case 'MPScanStaticSet': return MPScanStaticSet.from_prev_calc(
                prev_calc_dir=dir_path, 
                user_incar_settings=incar_corrections, 
                user_kpoints_settings=kpoints_corrections, 
                user_potcar_settings=potcar_corrections, 
                user_potcar_functional=potcar_functional_correction
            )
            case str(): raise ValueError(f'Provided string is not a valid preset name ({preset}).')
            case _: raise TypeError(f'"preset" arg expected a str type, got {type(preset)} instead.')

    match preset:
        case 'MPStaticSet': return MPStaticSet(
            structure=struct_or_path, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MatPESStaticSet': return MatPESStaticSet(
            structure=struct_or_path, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case 'MPScanStaticSet': return MPScanStaticSet(
            structure=struct_or_path, 
            user_incar_settings=incar_corrections, 
            user_kpoints_settings=kpoints_corrections, 
            user_potcar_settings=potcar_corrections, 
            user_potcar_functional=potcar_functional_correction
        )
        case str(): raise ValueError(f'Provided string is not a valid preset name ({preset}).')
        case _: raise TypeError(f'"preset" arg expected a str type, got {type(preset)} instead.')

########################################
def vasp_relaxation_settings(
        structure: SiteCollection, 
        preset: PMGRelaxSetType = 'MITRelaxSet', 
        user_corrections: Optional[Dict] = None
    ) -> VaspInput:
    '''
    Setup VASP inputs for a given structure using one of the pymatgen relaxation presets.

    Parameters:
        structure (SiteCollection):     The structure to write VASP inputs for.

        preset (str):                   The pymatgen preset to use for VASP inputs initialization.

        user_corrections (dict):        User defined settings. It allows to override some of the preset 
                                        INCAR, KPOINTS or POTCAR settings if necessary. Defaults to None.
    '''

    assert isinstance(structure, SiteCollection), \
    '"structure" argument format not supported. \
    It must be an instance of the SiteCollection class or one of its subclasses.'
    
    assert preset in PMGRelaxSet, \
    f'"preset" argument not recognized. \
    It must be one of the allowed pymatgen relaxation presets:\n\
    {PMGRelaxSet}.'

    assert isinstance(user_corrections, dict) or user_corrections is None, \
    'user_incar_settings must be a dict or None'

    vasp_input = _RelaxSet_init(
        structure=structure, 
        preset=preset, 
        corrections=user_corrections
    ).get_vasp_input()

    return vasp_input

########################################

def vasp_static_settings(
        structure: Optional[SiteCollection] = None, 
        preset: PMGStaticSetType = 'MPStaticSet', 
        from_prev_calc: bool = False, 
        prev_calc_dir: Optional[PathLike] = None, 
        user_corrections: Optional[dict] = None
    ) -> VaspInput:
    '''
    Setup VASP inputs for a given structure using one of the pymatgen static presets.

    Parameters:
        structure (SiteCollection):     The structure to write VASP inputs for.

        preset (str):                   The pymatgen preset to use for VASP inputs initialization.

        from_prev_calc (bool):          Whether to get final structure, INCAR and KPOINTS settings 
                                        from a previous VASP run.
                                        INCAR tags will still be managed to fit a static calculation 
                                        if previous run is a relaxation.
                                        For the sake of consistency, it is recommended to use the 
                                        static preset corresponding to previous relaxation preset in 
                                        this case (e.g. MPStaticSet for a relaxation with MPRelaxSet).
                                        If set to True, a directory to extract data from must be provided, 
                                        and structure argument is ignored. Defaults to False.

        prev_calc_dir (str|Path):       Directory to extract previous VASP run data from when from_prev_calc is True.
                                        If from_prev_calc is False, this argument is ignored.

        user_corrections (dict):        User defined settings. It allows to override some of the preset INCAR, KPOINTS 
                                        or POTCAR settings if necessary. Defaults to None.
    '''

    assert isinstance(structure, SiteCollection) or structure is None, \
    '"structure" argument format not supported. \
    It must be an instance of the SiteCollection class or one of its subclasses.'
    
    assert preset in PMGStaticSet, \
    f'"preset" argument not recognized. \
    It must be one of the allowed pymatgen static presets:\n\
    {PMGStaticSet}.'

    assert isinstance(user_corrections, dict) or user_corrections is None, \
    'user_corrections must be a dict or None'

    if not from_prev_calc:
        vasp_input = _StaticSet_init(
            struct_or_path=structure, 
            preset=preset, 
            corrections=user_corrections
        ).get_vasp_input()
    
    else:
        assert isinstance(prev_calc_dir, PathLike), \
        'from_prev_calc was set to True, prev_calc_dir must be provided as str or Path object.'

        prev_calc_dir = Path(prev_calc_dir)

        assert prev_calc_dir.is_dir(), \
        f'Prev_calc_dir: {prev_calc_dir} is not a valid directory.'

        vasp_input = _StaticSet_init(
            struct_or_path=prev_calc_dir, 
            from_prev_calc=from_prev_calc, 
            preset=preset, 
            corrections=user_corrections
        ).get_vasp_input()

    return vasp_input

########################################

def extract_vasp_data_for_convex_hull(
        struct_dir: PathLike = '.', 
        ignore_file: str = 'rejected.txt'
) -> Tuple[str, Dict[str, Union[Structure, float]]]:
    '''
    Extracts VASP data from a previous run for one structure. 
    Keeps only relevant data for relative stability calculation.

    Parameters:
        struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.

        ignore_file (str):      Checks whether the provided file name exists in structure directory.
                                Structure directories containing this file will return None.
                                This parameter permits the filtration of structures that did not pass
                                previous screening steps.
        
    Returns:
        Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict, 
                                containing following data, used in stability calculation:
                                    - structure chemical composition as Composition object, 
                                    - structure final energy (in eV).
    '''

    assert isinstance(struct_dir, PathLike)
    assert Path(struct_dir).is_dir()

    struct_dir: Path = Path(struct_dir)
    files            = set(file.name for file in struct_dir.iterdir())

    if ignore_file in files: return {}

    struct_name  = struct_dir.name
    contcar_path = Path(struct_dir / 'CONTCAR')
    oszicar_path = Path(struct_dir / 'OSZICAR')

    structure    = Poscar.from_file(contcar_path).structure
    composition  = Composition(structure.formula)
    final_energy_eV: float = Oszicar(oszicar_path).final_energy
    struct_dict = {
        'composition': composition, 
        'final_energy': final_energy_eV
    }
    struct_data = (struct_name, struct_dict)

    return struct_data

########################################

def extract_vasp_data_for_delta_sol(
        struct_dir: PathLike = '.', 
        ignore_file: str = 'rejected.txt'
) -> Tuple[str, Dict[str, Union[Structure, Chgcar, float]]]:
    '''
    Extracts VASP data from a previous run for one structure. 
    Keeps only relevant data for Δ-Sol method.

    Parameters:
        struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.

        ignore_file (str):      Checks whether the provided file name exists in structure directory.
                                Structure directories containing this file will return None.
                                This parameter permits the filtration of structures that did not pass
                                previous screening steps.
        
    Returns:
        Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict, 
                                containing following data, used in Δ-Sol method:
                                    - structure itself, 
                                    - its CHGCAR file (to modify charge density), 
                                    - its final energy (in eV), used as E(N0).
    '''

    assert isinstance(struct_dir, PathLike)
    assert Path(struct_dir).is_dir()

    struct_dir: Path = Path(struct_dir)
    files            = set(file.name for file in struct_dir.iterdir())

    if ignore_file in files: return {}

    struct_name  = struct_dir.name
    contcar_path = Path(struct_dir / 'CONTCAR')
    chgcar_path  = Path(struct_dir / 'CHGCAR')
    oszicar_path = Path(struct_dir / 'OSZICAR')

    structure       = Poscar.from_file(contcar_path).structure
    chgcar          = Chgcar.from_file(chgcar_path)
    final_energy_eV = Oszicar(oszicar_path).final_energy

    struct_dict = {
        'structure': structure, 
        'CHGCAR': chgcar, 
        'final_energy': final_energy_eV
    }
    struct_data = (struct_name, struct_dict)

    return struct_data

########################################

def batch_extract_vasp_data(
        method: Literal['convex_hull', 'delta_sol'], 
        base_dir: PathLike = '.', 
        ignore_file: str = 'rejected.txt', 
        workers: int = 1
) -> Dict[str, Dict[str, Any]]:
    '''
    Extracts VASP data from a previous run for each structure directory in given directory.
    Keeps only data that are useful according to given method.

    Parameters:
        method (str):           Name of the method that will use the data, used to know 
                                which data should be extracted.
                                Actual methods supported: 'convex_hull', 'delta_sol'.

        base_dir (str|Path):    Directory containing structures subdirs to extract data from.

        ignore_file (str):      Checks whether the provided file name exists in each 
                                subdirectory. Structure directories containing this file 
                                will not be taken into account.
                                This parameter permits the filtration of structures that did 
                                not pass previous steps.

        workers (int):          Number of parallel processes to spawn.

    Returns:
        Dict[str, Dict]:        Dict with structure directory names as keys, 
                                and a dict containing useful data according to chosen method 
                                for corresponding structure as values.

                                Data returned for 'convex_hull' method:
                                    - composition of the formula unit, 
                                    - final energy of the relaxation in eV.
                                
                                Data returned for 'delta_sol' method:
                                    - structure itself, 
                                    - its CHGCAR file (to modify charge density), 
                                    - final energy of the relaxation in eV, used as E(N0).
    '''

    assert isinstance(base_dir, PathLike)
    assert Path(base_dir).is_dir()
    assert isinstance(ignore_file, str)
    assert isinstance(workers, int) and workers >= 1

    def is_struct_dir(path: Path) -> bool:
        return path.is_dir() and path.name[0].isdecimal()
    
    match method:
        case "convex_hull":
            set_vasp_extractor = partial(extract_vasp_data_for_convex_hull, ignore_file=ignore_file)
        case "delta_sol":
            set_vasp_extractor = partial(extract_vasp_data_for_delta_sol, ignore_file=ignore_file)
        case _: raise NotImplementedError(f"Provided method ({method}) is not supported.")

    base_dir           = Path(base_dir)
    structs_dir_list   = list(filter(is_struct_dir, base_dir.iterdir()))
    nbr_structs        = len(structs_dir_list)
    chunksize          = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)

    structs_data_list  = list(filter(
        None, 
        process_map(
            set_vasp_extractor, 
            structs_dir_list, 
            max_workers=workers, 
            chunksize=chunksize, 
            desc='Extracting infos from previous VASP output'
        )
    ))
    structs_data = dict(structs_data_list)

    return structs_data

########################################

def chgcar_density_switch(chgcar: Chgcar, delta: float):
    '''
    Modifies provided Chgcar object's charge density by +/- delta, 
    then returns the two resulting Chgcar objects.
    '''

    e = elementary_charge # Exact value in Coulomb
    chgcar_plus, chgcar_minus = chgcar.copy(), chgcar.copy()
    chgcar_plus.data['total'] += delta*e
    chgcar_minus.data['total'] -= delta*e

    return chgcar_plus, chgcar_minus

########################################

def struct_charge_switch(structure: Structure, new_charge: float):
    '''
    Modifies overall charge of provided structure by +/- charge, 
    then returns the two resulting structures.

    Parameters:
        structure (Structure):  Neutral base structure on which charges will be added.

        new_charge (float):     Value of the charge to apply.
    
    Returns:
        Two copies of the input structure, with a positive and negative charge, respectively.
    '''

    pos_struct, neg_struct = structure.copy(), structure.copy()
    pos_struct.set_charge(structure.charge + new_charge)
    neg_struct.set_charge(structure.charge - new_charge)

    return pos_struct, neg_struct

########################################

def delta_sol_inputs_init(
        structs_data: dict, 
        preset: PMGStaticSetType = 'MPStaticSet', 
        user_corrections: Optional[Dict] = None, 
        with_uncertainties: bool = False
    ) -> Tuple[List]:
    '''
    Initialize Vasp static input sets for structures with N0 - n and N0 + n electrons per cell, 
    used in band gap calculations with Δ-Sol method.

    Reference:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (reference 32 in screening_pipeline/Bibliography, values in Table I)

    Parameters:
        structs_data (dict):        Dict containing relevant VASP data from a previous relaxation, 
                                    as provided by the batch_extract_vasp_data function.
        
        preset (PMGStaticSet):      One of the pymatgen static VASP preset labeled 'StaticSet'.

        with_uncertainties (bool):  Whether to also compute uncertainty boundaries of Δ-Sol method.
                                    Defaults to False.
    
    Returns:
        Tuple[List[str], List[VaspInput]]:  
        A list of the subdirectory names corresponding to Δ-Sol runs, 
        and a list of the VaspInput objects for these runs, both in same order.
    '''
    
    def _inputs_init(struct_tuple: Tuple[str, Dict]) -> List[Tuple[str, VaspInput]]:

        # Initialize input sets
        name      = struct_tuple[0]
        data      = struct_tuple[1]
        structure = data['structure']
        N_val     = get_all_valence_electrons(structure)
        run_set   = vasp_static_settings(structure, preset=preset, user_corrections=user_corrections)

        # Search for the right N* parameter to use
        pot_func = run_set.get('POTCAR_FUNCTIONAL', 'PBE')

        if 'LDA' in pot_func: delta_sol_functional = 'LDA'
        elif 'PBE' in pot_func: delta_sol_functional = 'PBE'
        elif 'AM05' in pot_func: delta_sol_functional = 'AM05'
        else:
            raise NotImplementedError(
                "Provided POTCAR functional is not implemented for Δ-Sol method. \
                Recognized functionals are the ones having 'LDA', 'PBE' or 'AM05' in the name."
            )

        struct_runs_list = []
        n_star_types     = ('BEST', 'MIN', 'MAX') if with_uncertainties else ('BEST',)

        # Compute relevant n = N_val / N* 
        for n_star_type in n_star_types:

            n_ratio = get_delta_sol_el_ratio(
                structure=structure, 
                dft_functional=delta_sol_functional, 
                n_star_type=n_star_type
            )

            run_plus       = run_set.copy().update({'NELECT': N_val + n_ratio})
            run_minus      = run_set.copy().update({'NELECT': N_val - n_ratio})
            run_plus_path  = Path('_'.join((name , n_star_type.lower(), 'plus')))
            run_minus_path = Path('_'.join((name , n_star_type.lower(), 'minus')))

            struct_runs_list += [(run_plus_path, run_plus), (run_minus_path, run_minus)]

        return struct_runs_list

    input_data   = flatten(list(starmap(_inputs_init, structs_data.items())))
    subdirs_list = [tup[0] for tup in input_data]
    inputs_list  = [tup[1] for tup in input_data]

    return subdirs_list, inputs_list

########################################
