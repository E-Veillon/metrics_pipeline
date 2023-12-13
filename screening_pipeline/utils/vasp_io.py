'''
Functions to write VASP input files, launch VASP calculations and manage VASP output files.
'''


########################################
# SYSTEM I/O MODULES

from typing import Optional, Dict, List, Union, Sequence, Tuple, Literal
from pathlib import Path

########################################
# OPTIMIZATION MODULES

from scipy.constants import elementary_charge
from itertools import cycle, chain, starmap
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure, SiteCollection
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
#TODO: This function might be useless, delete it if no use at the end of pipeline devpt
def vasp_input_files_settings(
        structure: Structure, 
        /, *, 
        use_mit_set: bool = True, 
        corrected: bool = True, 
        config_dict: Optional[dict] = None, 
        user_incar_settings: Optional[dict] = None, 
        user_kpoints_settings: Optional[dict] = None, 
        user_potcar_settings: Optional[dict] = None
    ) -> VaspInput:
    '''
    Builds a VaspInput object, with a standard preset option corresponding to the MIT high throughput material screening project.
        
    Reference of the preset:
        A. Jain, G. Hautier, C.J. Moore, S.P. Ong,
        C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder,
        Computational Materials Science, 50, 2295-2310 (2011)
        (reference 14 in screening_pipeline/Bibliography/)

    Parameters:
        structure (Structure): The structure to write VASP input files for.

        config_dict (dict):             VASP input parameters provided as a dict.
                                        Input file names are the keys and dicts containing file specific parameters are the values.

        use_mit_set (bool):             Whether to use the MIT high throughput project's preset or not.
                                        If True, config_dict should be set to None, but input parameters can be ajusted with appropriate modifier arguments (default).
                                        If False, you need to provide a full config_dict.

        user_incar_settings (dict):     User INCAR settings. It allows to override some of the standard INCAR tags if necessary.
                                        Defaults to None.

        user_kpoints_settings (dict):   User KPOINTS settings. It allows to override the standard Kpoints mesh if necessary.
                                        Defaults to None.

        user_potcar_settings (dict):    User POTCAR settings. It allows to override the standard POTCAR settings, although it is not recommended.
                                        Defaults to None.
    Returns:
        A DictSet object that uses the write_input method to write set input files in given directory, ready for calculation.
    '''

    assert isinstance(structure, Structure), 'Provided structure format is not supported. Please provide a PyMatGen Structure object or one of its subclasses'
    
    if not use_mit_set:

        assert isinstance(config_dict, dict), 'use_mit_set was set to False, a full config_dict has to be provided.'
        assert 'INCAR' in config_dict.keys(), 'No INCAR configuration detected in provided config_dict.'
        assert ('KPOINTS' in config_dict.keys()) or ('KSPACING' in config_dict['INCAR']), 'No k-points configuration detected, either by a KPOINTS key in config_dict or by a KSPACING INCAR tag.'
        assert 'POTCAR' in config_dict.keys(), 'No POTCAR configuration detected in provided config_dict.'

        return DictSet(structure, config_dict).get_vasp_input()

    MITRelaxSet_corrections_dict = {}

    if corrected:
        MITRelaxSet_corrections_dict = {'INCAR': _MITRelaxSet_INCAR_corrections(structure.num_sites)}

    def has_only_string_keys(dict: Dict) -> bool:
        return all(map(isinstance,dict.keys(),cycle((str,))))

    if user_incar_settings is not None:
        assert isinstance(user_incar_settings, dict), 'INCAR modifications should be provided as a dict.'
        assert has_only_string_keys(user_incar_settings), 'All provided INCAR tags should be strings.'

        for incar_tag, tag_value in user_incar_settings:
            MITRelaxSet_corrections_dict['INCAR'][incar_tag] = tag_value

    if user_kpoints_settings is not None:
        assert isinstance(user_kpoints_settings, dict), 'KPOINTS modifications should be provided as a dict.'
        assert has_only_string_keys(user_kpoints_settings), 'All provided KPOINTS modifications keys should be strings.'

        MITRelaxSet_corrections_dict['KPOINTS'] = user_kpoints_settings

    if user_potcar_settings is not None:
        assert isinstance(user_potcar_settings, dict), 'POTCAR modifications should be provided as a dict.'
        assert has_only_string_keys(user_potcar_settings), 'All provided POTCAR modifications keys should be strings.'

        MITRelaxSet_corrections_dict['POTCAR'] = user_potcar_settings

    vasp_input = MITRelaxSet(
        structure, 
        user_incar_settings=MITRelaxSet_corrections_dict['INCAR'], 
        user_kpoints_settings=MITRelaxSet_corrections_dict['KPOINTS'], 
        user_potcar_settings=MITRelaxSet_corrections_dict['POTCAR']
    ).get_vasp_input()

    return vasp_input

########################################
#TODO: This function might be useless, delete it if no use at the end of pipeline devpt
def batch_write_MITRelaxSet_inputs(
        structures: List[Structure], 
        corrected: bool = True, 
        workers: int = 1, 
        **kwargs
    ) -> List[VaspInput]:
    '''
    Creates VaspInput objects for several structures.
    Uses the MITRelaxSet preset from pymatgen.

    Reference of the preset:
        A. Jain, G. Hautier, C.J. Moore, S.P. Ong,
        C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder,
        Computational Materials Science, 50, 2295-2310 (2011)
        (reference 14 in screening_pipeline/Bibliography/)

    Parameters:
        structures (List[Structure]):   The structures to write VASP inputs for.
        
        corrected (bool):               The preset implemented in pymatgen uses slightly different parameters
                                        compared to the ones cited in the reference article. Setting this to True
                                        corrects the preset to better fit the reference. Defaults to True.
        
        kwargs:                         Any user settings supported by MITRelaxSet.
    '''
    
    assert all(map(isinstance, structures, cycle((Structure,)))), \
    'Some of the structures provied are not Structure objects'

    user_incar_settings   = kwargs.get('user_incar_settings', None)
    user_kpoints_settings = kwargs.get('user_kpoints_settings', None)
    user_potcar_settings  = kwargs.get('user_potcar_settings', None)

    assert isinstance(user_incar_settings, (dict, None)), \
    'user_incar_settings must be a dict or None'
    assert isinstance(user_kpoints_settings, (dict, None)), \
    'user_kpoints_settings must be a dict or None'
    assert isinstance(user_potcar_settings, (dict, None)), \
    'user_potcar_settings must be a dict or None'

    vasp_input_files_settings_part = partial(
        vasp_input_files_settings, 
        use_mit_set = True, 
        corrected=corrected, 
        user_incar_settings=user_incar_settings, 
        user_kpoints_settings=user_kpoints_settings, 
        user_potcar_settings=user_potcar_settings
    )
    
    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    vasp_inputs = list(
        process_map(
            vasp_input_files_settings_part,
            structures, 
            workers=workers, 
            chunksize=chunksize, 
            desc='Converting structures to VASP inputs'
        )
    )
    return vasp_inputs

########################################

def vasp_launcher(vasp_input: VaspInput, path: PathLike) -> None:

    assert isinstance(vasp_input, VaspInput)
    assert isinstance(path, PathLike)
    path     = Path(path)
    calc_dir = path
    out_file = path / "vasp.out"
    err_file = path / "vasp.err"
    vasp_input.run_vasp(run_dir=calc_dir, output_file=out_file, err_file=err_file)

########################################

def vasp_batch_launch(
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
        vasp_inputs ([VaspInput]):  The objects defining how to write VASP input files in each subdirectory.

        base_dir (str|Path):        The directory where the subdirs should be created.

        subdir_names ([str|Path]):  The names or subpaths for created subdirectories.
                                    Note that its length must match the length of vasp_inputs.

        workers (int):              The number of parallel processes to spawn.
                                    Defaults to 1.
    '''

    assert all(isinstance(vasp_input, VaspInput) for vasp_input in vasp_inputs)
    assert isinstance(base_dir, PathLike)
    
    base_dir = Path(base_dir)

    assert base_dir.is_dir()
    assert all(isinstance(subdir_name, PathLike) for subdir_name in subdir_names)
    assert len(vasp_inputs) == len(subdir_names)
    assert isinstance(workers, int) and workers >= 1

    subpaths_list = batch_add_new_dirs(base_dir=base_dir, new_subdirs=subdir_names)
    inputs_list   = list(zip(vasp_inputs, subpaths_list))

    process_map(
        vasp_launcher, 
        inputs_list, 
        workers=workers, 
        chunksize=1
    )

########################################

def vasp_relaxation_settings(
        structure: SiteCollection, 
        preset: PMGRelaxSet = 'MITRelaxSet', 
        user_incar_settings: Optional[dict] = None, 
        user_kpoints_settings: Optional[dict] = None, 
        user_potcar_settings: Optional[dict] = None
    ) -> VaspInput:
    '''
    Setup VASP inputs for a given structure using one of the pymatgen relaxation presets.

    Parameters:
        structure (SiteCollection):     The structure to write VASP inputs for.

        preset (str):                   The pymatgen preset to use for VASP inputs initialization.

        user_incar_settings (dict):     User INCAR settings. It allows to override some of the standard INCAR tags if necessary.
                                        Defaults to None.

        user_kpoints_settings (dict):   User KPOINTS settings. It allows to override the standard Kpoints setup if necessary.
                                        Defaults to None.

        user_potcar_settings (dict):    User POTCAR settings. It allows to override the standard POTCAR settings, although it is not recommended.
                                        Defaults to None.
    '''
    
    from pymatgen.io.vasp.sets import   MITRelaxSet, MPRelaxSet, MPScanRelaxSet, \
                                        MPHSERelaxSet, MPMetalRelaxSet, MVLRelax52Set, \
                                        MVLScanRelaxSet
    
    allowed_presets = {
        'MITRelaxSet': MITRelaxSet, 
        'MPRelaxSet': MPRelaxSet, 
        'MPScanRelaxSet': MPScanRelaxSet, 
        'MPHSERelaxSet': MPHSERelaxSet, 
        'MPMetalRelaxSet': MPMetalRelaxSet, 
        'MVLRelax52Set': MVLRelax52Set, 
        'MVLScanRelaxSet': MVLScanRelaxSet
    }

    assert isinstance(structure, SiteCollection), '''
    "structure" argument format not supported.
    It must be an instance of the SiteCollection class or one of its subclasses.'''
    
    assert preset in allowed_presets.keys(), '''
    "preset" argument not recognized.
    It must be one of the allowed pymatgen relaxation presets.'''

    assert isinstance(user_incar_settings, (dict, None)), \
    'user_incar_settings must be a dict or None'

    assert isinstance(user_kpoints_settings, (dict, None)), \
    'user_kpoints_settings must be a dict or None'

    assert isinstance(user_potcar_settings, (dict, None)), \
    'user_potcar_settings must be a dict or None'

    if preset == 'MITRelaxSet':
        MIT_INCAR_corrections = _MITRelaxSet_INCAR_corrections(structure.num_sites)

        if user_incar_settings is not None:
            MIT_INCAR_corrections.update(user_incar_settings)
            user_incar_settings.update(MIT_INCAR_corrections)

    preset_obj: DictSet = allowed_presets[preset]

    vasp_input = preset_obj(
        structure=structure, 
        user_incar_settings=user_incar_settings, 
        user_kpoints_settings=user_kpoints_settings, 
        user_potcar_settings=user_potcar_settings
    ).get_vasp_input()

    return vasp_input

########################################

def vasp_static_settings(
        structure: SiteCollection, 
        preset: PMGStaticSet = 'MPStaticSet', 
        from_prev_calc: bool = False, 
        prev_calc_dir: Optional[PathLike] = None, 
        user_incar_settings: Optional[dict] = None, 
        user_kpoints_settings: Optional[dict] = None, 
        user_potcar_settings: Optional[dict] = None
    ) -> VaspInput:
    '''
    Setup VASP inputs for a given structure using one of the pymatgen static presets.

    Parameters:
        structure (SiteCollection):     The structure to write VASP inputs for.

        preset (str):                   The pymatgen preset to use for VASP inputs initialization.

        from_prev_calc (bool):          Whether to get final structure, INCAR and KPOINTS settings from a previous VASP run.
                                        INCAR tags will still be managed to fit a static calculation if previous run is a relaxation.
                                        For the sake of consistency, it is recommended to use the static preset corresponding to
                                        previous relaxation preset in this case (e.g. MPStaticSet for a relaxation with MPRelaxSet).
                                        If set to True, a directory to extract data from must be provided. Defaults to False.

        prev_calc_dir (str|Path):       Directory to extract previous VASP run data from when from_prev_calc is True.
                                        If from_prev_calc is False, this argument is ignored.

        user_incar_settings (dict):     User INCAR settings. It allows to override some of the standard INCAR tags if necessary.
                                        Defaults to None.

        user_kpoints_settings (dict):   User KPOINTS settings. It allows to override the standard Kpoints setup if necessary.
                                        Defaults to None.

        user_potcar_settings (dict):    User POTCAR settings. It allows to override the standard POTCAR settings, although it is not recommended.
                                        Defaults to None.
    '''

    from pymatgen.io.vasp.sets import MPStaticSet, MatPESStaticSet, MPScanStaticSet

    allowed_presets = {
        'MPStaticSet': MPStaticSet, 
        'MatPESStaticSet': MatPESStaticSet, 
        'MPScanStaticSet': MPScanStaticSet
    }

    assert isinstance(structure, SiteCollection), '''
    "structure" argument format not supported.
    It must be an instance of the SiteCollection class or one of its subclasses.'''
    
    assert preset in allowed_presets.keys(), '''
    "preset" argument not recognized.
    It must be one of the allowed pymatgen static presets.'''

    assert isinstance(user_incar_settings, (dict, None)), \
    'user_incar_settings must be a dict or None'

    assert isinstance(user_kpoints_settings, (dict, None)), \
    'user_kpoints_settings must be a dict or None'

    assert isinstance(user_potcar_settings, (dict, None)), \
    'user_potcar_settings must be a dict or None'

    preset_obj: DictSet = allowed_presets[preset]

    if not from_prev_calc:
        vasp_input = preset_obj(
            structure=structure, 
            user_incar_settings=user_incar_settings, 
            user_kpoints_settings=user_kpoints_settings, 
            user_potcar_settings=user_potcar_settings
        ).get_vasp_input()
    
    else:
        assert isinstance(prev_calc_dir, PathLike), \
        'from_prev_calc was set to True, prev_calc_dir must be provided as str or Path object'

        prev_calc_dir = Path(prev_calc_dir)

        assert prev_calc_dir.is_dir(), \
        'Provided prev_calc_dir is not a valid directory'

        vasp_input = preset_obj.from_prev_calc(
            prev_calc_dir=prev_calc_dir, 
            user_incar_settings=user_incar_settings, 
            user_kpoints_settings=user_kpoints_settings, 
            user_potcar_settings=user_potcar_settings
        ).get_vasp_input()

    return vasp_input

########################################

def extract_vasp_data_from_prev_calc(
        struct_dir: PathLike = '.', 
        ignore_file: str = 'rejected.txt'
) -> Tuple[str, Dict[Structure, Chgcar, float]]:
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
                                containing following useful data:
                                    - structure itself, 
                                    - its CHGCAR file (to modify charge density), 
                                    - its final energy (in eV), used as E(N0).
    '''
    assert isinstance(struct_dir, PathLike)

    struct_dir: Path = Path(struct_dir)

    assert struct_dir.is_dir()

    files = set(file.name for file in struct_dir.iterdir())

    if ignore_file in files:
        return None

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
        base_dir: PathLike = '.', 
        ignore_file: str = 'rejected.txt', 
        workers: int = 1
) -> Dict[str, Dict[str, Union[Structure, Chgcar, float]]]:
    '''
    Extracts VASP data from a previous run for each structure directory in given directory.
    Keeps only data that are useful for Δ-Sol method.

    Parameters:
        base_dir (str|Path):    Directory containing structures subdirs to extract data from.

        ignore_file (str):      Checks whether the provided file name exists in each subdirectory.
                                Structure directories containing this file will not be taken into account.
                                This parameter permits the filtration of structures that did not pass
                                previous steps.

        workers (int):          Number of parallel processes to spawn.

    Returns:
        Dict[str: Dict]:        Dict with structure directory names as keys, 
                                and a dict containing following data for corresponding structure as values:
                                    - structure itself, 
                                    - its CHGCAR file (to modify charge density), 
                                    - its final energy (in eV), used as E(N0).
    '''

    assert isinstance(base_dir, PathLike)

    base_dir: Path = Path(base_dir)

    assert base_dir.is_dir()

    def is_directory(path: Path) -> bool:
        return path.is_dir()

    set_vasp_extractor = partial(extract_vasp_data_from_prev_calc, ignore_file=ignore_file)
    structs_dir_list   = list(filter(is_directory, base_dir.iterdir()))
    nbr_structs        = len(structs_dir_list)
    chunksize          = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)

    structs_data_list  = list(filter(
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
            dft_functional: Literal['LDA', 'PBE', 'AM05'] = 'PBE', 
            n_star_type: Literal['MIN', 'BEST', 'MAX'] = 'BEST'
        ) -> List[Tuple[str, VaspInput]]:

        # Produce E(N0 + n) and E(N0 - n)'s CHGCAR files
        structure       = data['structure']
        chgcar          = data['CHGCAR']
        data['n_ratio'] = get_delta_sol_el_ratio(structure, dft_functional, n_star_type)
        data['CHGCAR_plus'], data['CHGCAR_minus'] = chgcar_density_switch(chgcar, data['n_ratio'])

        # Prepare E(N0 + n) input set
        run_plus = vasp_static_settings(structure, preset=preset)
        run_plus.update({'CHGCAR': data['CHGCAR_plus']})
        run_plus_path = Path('_'.join(name , 'plus'))

        # Prepare E(N0 - n) input set
        run_minus = vasp_static_settings(structure, preset=preset)
        run_minus.update({'CHGCAR': data['CHGCAR_minus']})
        run_minus_path = Path('_'.join(name , 'minus'))

        struct_list = [(run_plus_path, run_plus), (run_minus_path, run_minus)]

        return struct_list

    set_inputs_init = partial(
        inputs_init, 
        preset=preset, 
        dft_functional=dft_functional, 
        n_star_type=n_star_type
    )
    structs_tuples = list(chain(starmap(set_inputs_init, structs_data.items())))
    subdirs_list   = [tup[0] for tup in structs_tuples]
    inputs_list    = [tup[1] for tup in structs_tuples]

    return inputs_list, subdirs_list