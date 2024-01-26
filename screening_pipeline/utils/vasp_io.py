'''
Functions to write VASP input files, launch VASP calculations and manage VASP output files.
'''


########################################
# SYSTEM I/O MODULES

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
from pymatgen.io.vasp import VaspInput
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Chgcar, Oszicar
#TODO: The next import might be useless. Delete it if related functions are not used until end of pipeline devpt.
from pymatgen.io.vasp.sets import DictSet, MITRelaxSet

########################################
# LOCAL MODULES

from screening_pipeline.utils.fitted_values import U_VALUES
from screening_pipeline.utils.paths import batch_add_new_dirs
from screening_pipeline.utils.periodic_table import get_delta_sol_el_ratio

########################################
# TYPE ALIASES

PathLike = Union[str, Path]
PMGRelaxSet = Literal[
    'MITRelaxSet', 
    'MPRelaxSet', 
    'MPScanRelaxSet', 
    'MPHSERelaxSet', 
    'MPMetalRelaxSet', 
    'MVLRelax52Set', 
    'MVLScanRelaxSet'
]
PMGStaticSet = Literal[
    'MPStaticSet', 
    'MatPESStaticSet', 
    'MPScanStaticSet'
]

########################################
# LOCAL FUNCTIONS

def _MITRelaxSet_INCAR_corrections(number_of_sites: int) -> Dict:
    '''
    Corrects errors and imprecisions found in Pymatgen MITRelaxSet VASP preset's INCAR tags.
    For example, some tags should depend on the number of atoms per unit cell of the structure.

    parameters:
        number_of_sites (int): number of atoms in a unit cell.
    
    Returns:
        A dictionnary containing the INCAR tags corrections.
    '''

    corrected_EDIFF = float(5e-5)*number_of_sites
    corrected_ENCUT = 520 # To be modified according to ENMAX value (ENCUT = 1.3*ENMAX)
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
    vasp_input.run_vasp(run_dir=calc_dir, vasp_cmd=vasp_exe_list, output_file=out_file, err_file=err_file)

########################################

def _vasp_launcher_wrapper(args_tuple: Tuple[VaspInput, PathLike]):
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
        _vasp_launcher_wrapper, 
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
        incar_corrections = _MITRelaxSet_INCAR_corrections(structure.num_sites)
    
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
        preset: PMGRelaxSet = 'MITRelaxSet', 
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
    
    allowed_presets = {
        'MITRelaxSet', 
        'MPRelaxSet', 
        'MPScanRelaxSet', 
        'MPHSERelaxSet', 
        'MPMetalRelaxSet', 
        'MVLRelax52Set', 
        'MVLScanRelaxSet'
    }

    assert isinstance(structure, SiteCollection), '''
    "structure" argument format not supported.
    It must be an instance of the SiteCollection class or one of its subclasses.'''
    
    assert preset in allowed_presets, f'''
    "preset" argument not recognized.
    It must be one of the allowed pymatgen relaxation presets:
    {allowed_presets}'''

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
        structure: Optional[SiteCollection], 
        preset: PMGStaticSet = 'MPStaticSet', 
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

    allowed_presets = {
        'MPStaticSet', 
        'MatPESStaticSet', 
        'MPScanStaticSet'
    }

    assert isinstance(structure, SiteCollection) or structure is None, '''
    "structure" argument format not supported.
    It must be an instance of the SiteCollection class or one of its subclasses.'''
    
    assert preset in allowed_presets.keys(), f'''
    "preset" argument not recognized.
    It must be one of the allowed pymatgen static presets:
    {allowed_presets}'''

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
        'from_prev_calc was set to True, prev_calc_dir must be provided as str or Path object'

        prev_calc_dir = Path(prev_calc_dir)

        assert prev_calc_dir.is_dir(), \
        'Provided prev_calc_dir is not a valid directory'

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

    structure    = Poscar.from_file(contcar_path).structure
    chgcar       = Chgcar.from_file(chgcar_path)
    final_energy_eV: float = Oszicar(oszicar_path).final_energy
    #final_energy_eV_per_at = final_energy_eV / float(structure.num_sites)
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

    def is_directory(path: Path) -> bool:
        return path.is_dir()
    
    match method:
        case "convex_hull":
            set_vasp_extractor = partial(extract_vasp_data_for_convex_hull, ignore_file=ignore_file)
        case "delta_sol":
            set_vasp_extractor = partial(extract_vasp_data_for_delta_sol, ignore_file=ignore_file)
        case _: raise NotImplementedError(f"Provided method ({method}) is not supported.")

    base_dir           = Path(base_dir)
    structs_dir_list   = list(filter(is_directory, base_dir.iterdir()))
    nbr_structs        = len(structs_dir_list)
    chunksize          = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)

    structs_data_list  = list(filter(
        None, 
        process_map(
            set_vasp_extractor, 
            structs_dir_list, 
            workers=workers, 
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

def struct_charge_switch(structure: Structure, delta: float):
    '''
    Modifies overall charge of provided structure by +/- delta, 
    then returns the two resulting structures.
    '''

    cation_struct, anion_struct = structure.copy(), structure.copy()
    cation_struct.set_charge(structure.charge + delta)
    anion_struct.set_charge(structure.charge - delta)
    return cation_struct, anion_struct

########################################

def delta_sol_inputs_init(
        structs_data: dict, 
        preset: PMGStaticSet = 'MPStaticSet', 
        user_corrections: Optional[Dict] = None, 
        dft_functional: Literal['LDA', 'PBE', 'AM05'] = 'PBE', 
        n_star_type: Literal['MIN', 'BEST', 'MAX'] = 'BEST'
    ) -> Tuple[List]:
    '''
    Initialize Vasp static input sets for structures with N0 - n and N0 + n electrons per cell, 
    used in band gap calculations with Δ-Sol method.

    Reference:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (reference 32 in screening_pipeline/Bibliography, values in Table I)

    Parameters:
        structs_data (dict):    Dict containing relevant VASP data from a previous relaxation, 
                                as provided by the batch_extract_vasp_data function.
        
        preset (PMGStaticSet):  One of the pymatgen static VASP preset labeled 'StaticSet'.

        dft_functional ('LDA'|'PBE'|'AM05'):    Choose one of the functionals supported by Δ-Sol method.
                                                Used to initialize N*. Defaults to 'PBE'.

        n_star_type ('MIN'|'BEST'|'MAX'):       Which N* to initialize for given functional.
                                                'BEST' is used for band gap estimation, 
                                                while 'MIN' and 'MAX' are for uncertainty measurements on said gap.
                                                Defaults to 'BEST'.
    
    Returns:
        Tuple[List]:    A list of VaspInput objects corresponding to Δ-Sol runs, 
                        and a list of the subdirectory names for these runs, in same order.
    '''
    
    def inputs_init(
            name: str, 
            data: Dict, 
            preset: PMGStaticSet = 'MPStaticSet', 
            user_corrections: Optional[Dict] = None, 
            dft_functional: Literal['LDA', 'PBE', 'AM05'] = 'PBE', 
            n_star_type: Literal['MIN', 'BEST', 'MAX'] = 'BEST'
        ) -> List[Tuple[str, VaspInput]]:

        # Produce E(N0 + n) and E(N0 - n)'s CHGCAR files
        structure       = data['structure']
        chgcar          = data['CHGCAR']
        data['n_ratio'] = get_delta_sol_el_ratio(structure, dft_functional, n_star_type)
        data['CHGCAR_plus'], data['CHGCAR_minus'] = chgcar_density_switch(chgcar, data['n_ratio'])

        # Prepare E(N0 + n) input set
        run_plus = vasp_static_settings(structure, preset=preset, user_corrections=user_corrections)
        run_plus.update({'CHGCAR': data['CHGCAR_plus']})
        run_plus_path = Path('_'.join((name , 'plus')))

        # Prepare E(N0 - n) input set
        run_minus = vasp_static_settings(structure, preset=preset, user_corrections=user_corrections)
        run_minus.update({'CHGCAR': data['CHGCAR_minus']})
        run_minus_path = Path('_'.join((name , 'minus')))

        struct_list = [(run_plus_path, run_plus), (run_minus_path, run_minus)]

        return struct_list

    set_inputs_init = partial(
        inputs_init, 
        preset=preset, 
        user_corrections=user_corrections, 
        dft_functional=dft_functional, 
        n_star_type=n_star_type
    )
    structs_tuples = list(chain(starmap(set_inputs_init, structs_data.items())))
    subdirs_list   = [tup[0] for tup in structs_tuples]
    inputs_list    = [tup[1] for tup in structs_tuples]

    return inputs_list, subdirs_list

########################################
