'''
Functions to write VASP input files, launch VASP calculations and manage VASP output files.
'''


########################################
# SYSTEM I/O MODULES

from typing import Optional, Dict, List, Union, Sequence
from pathlib import Path

########################################
# OPTIMIZATION MODULES

from itertools import cycle
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure, SiteCollection
from pymatgen.io.vasp import VaspInput
from pymatgen.io.vasp.sets import DictSet, MITRelaxSet

########################################
# LOCAL MODULES

from screening_pipeline.utils.fitted_values import U_VALUES
from screening_pipeline.utils.paths import batch_add_new_dirs

PathLike = Union[str, Path]

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

def vasp_launcher(vasp_input: VaspInput, path: PathLike) -> None:

    assert isinstance(vasp_input, VaspInput)
    assert isinstance(path, PathLike)
    path     = Path(path)
    calc_dir = path
    out_file = path / "vasp.out"
    err_file = path / "vasp.err"
    vasp_input.run_vasp(run_dir=calc_dir, output_file=out_file, err_file=err_file)

def vasp_batch_launch(
        vasp_inputs: Sequence[VaspInput], 
        base_dir: PathLike, 
        subdir_names: Sequence[str], 
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

        subdir_names ([str]):       The names of created subdirectories.
                                    Note that its length must match the length of vasp_inputs.

        workers (int):              The number of parallel processes to spawn.
                                    Defaults to 1.
    '''

    assert all(isinstance(vasp_input, VaspInput) for vasp_input in vasp_inputs)
    assert isinstance(base_dir, PathLike)
    
    base_dir = Path(base_dir)

    assert base_dir.is_dir()
    assert all(isinstance(subdir_name, str) for subdir_name in subdir_names)
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

def vasp_relaxation_settings(
        structure: SiteCollection, 
        preset: str = 'MITRelaxSet', 
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

def vasp_static_settings(
        structure: SiteCollection, 
        preset: str = 'MPStaticSet', 
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