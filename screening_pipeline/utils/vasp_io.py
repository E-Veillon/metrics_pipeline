'''
Functions to write VASP input files, launch VASP calculations and manage VASP output files.
'''


########################################
# TYPE HINTING

from typing import Optional, Dict

########################################
# OPTIMIZATION MODULES

from itertools import cycle

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import VaspInput
from pymatgen.io.vasp.sets import DictSet, MITRelaxSet

########################################
# LOCAL MODULES

from screening_pipeline.utils import U_VALUES


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

def vasp_input_files_settings(
        structure: Structure, 
        /, *, 
        use_mit_set: bool = True, 
        config_dict: Optional[dict] = None, 
        modified_incar: Optional[dict] = None, 
        modified_kpoints: Optional[dict] = None, 
        modified_potcar: Optional[dict] = None
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

            config_dict (dict):         VASP input parameters provided as a dictionnary with input file names as keys and dicts containing file specific parameters as values.

            use_mit_set (bool):         Whether to use the MIT high throughput project's preset or not.
                                        If True, config_dict should be set to None, but input parameters can be ajusted with appropriate modifier arguments (default).
                                        If False, you need to provide a full config_dict.

            modified_incar (dict):      User INCAR settings. It allows to override some of the standard INCAR tags if necessary.
                                        Defaults to None.

            modifieed_kpoints (dict):   User KPOINTS settings. It allows to override the standard Kpoints mesh if necessary.
                                        Defaults to None.

            modified_potcar (dict):     User POTCAR settings. It allows to override the standard POTCAR settings, although it is not recommended.
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

    MITRelaxSet_corrections_dict = {'INCAR': _MITRelaxSet_INCAR_corrections(structure.num_sites)}

    def has_only_string_keys(dict: Dict) -> bool:
        return all(map(isinstance,dict.keys(),cycle((str,))))

    if modified_incar is not None:
        assert isinstance(modified_incar, dict), 'INCAR modifications should be provided as a dict.'
        assert has_only_string_keys(modified_incar), 'All provided INCAR tags should be strings.'

        for incar_tag, tag_value in modified_incar:
            MITRelaxSet_corrections_dict['INCAR'][incar_tag] = tag_value

    if modified_kpoints is not None:
        assert isinstance(modified_kpoints, dict), 'KPOINTS modifications should be provided as a dict.'
        assert has_only_string_keys(modified_kpoints), 'All provided KPOINTS modifications keys should be strings.'

        for key, value in modified_kpoints:
            MITRelaxSet_corrections_dict['KPOINTS'][key] = value

    if modified_potcar is not None:
        assert isinstance(modified_potcar, dict), 'POTCAR modifications should be provided as a dict.'
        assert has_only_string_keys(modified_potcar), 'All provided POTCAR modifications keys should be strings.'

        for key, value in modified_potcar:
            MITRelaxSet_corrections_dict['POTCAR'][key] = value

    return MITRelaxSet(
        structure, 
        user_incar_settings=MITRelaxSet_corrections_dict['INCAR'], 
        user_kpoints_settings=MITRelaxSet_corrections_dict['KPOINTS'], 
        user_potcar_settings=MITRelaxSet_corrections_dict['POTCAR']
        ).get_vasp_input()

def vasp_launcher(vasp_input: VaspInput, path: str):
    calc_dir = path
    out_file = path + "vasp.out"
    err_file = path + "vasp.err"
    vasp_input.run_vasp(run_dir=calc_dir, output_file=out_file, err_file=err_file)