#!/usr/bin/python

from typing import Optional, Dict
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import DictSet, MITRelaxSet 

def _MITRelaxSet_INCAR_corrections(number_of_sites: int):
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
    corrected_LDAUU = {
        'F': {
            'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4.0, #'Cu': 4 -> 4.0
            'Fe': 4.0, 'Mn': 3.9, 'Mo': 3.5, 'Nb': 1.5, #'Mo': 4.38 -> 3.5 (according to the reference)
            'Ni': 6.0, 'Re': 2.0, 'Ta': 2.0, 'V': 3.1,  #'Ni': 6 -> 6.0, 'Re': 2 -> 2.0, 'Ta': 2 -> 2.0
            'W': 4.0
        }, 
        'O': {
            'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4.0, #'Cu': 4 -> 4.0
            'Fe': 4.0, 'Mn': 3.9, 'Mo': 3.5, 'Nb': 1.5, #'Mo': 4.38 -> 3.5
            'Ni': 6.0, 'Re': 2.0, 'Ta': 2.0, 'V': 3.1,  #'Ni': 6 -> 6.0, 'Re': 2 -> 2.0, 'Ta': 2 -> 2.0
            'W': 4.0                          
        }, 
        'S': {
            'Fe': 1.9, 'Mn': 2.5
        }}
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
        ):
    '''
    Builds a VASP input files generator object, with a standard preset possibility corresponding to the MIT high throughput material screening project.
        
        Reference of the preset:
            A. Jain, G. Hautier, C.J. Moore, S.P. Ong,
            C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder,
            Computational Materials Science, 50, 2295-2310 (2011)

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

        return DictSet(structure, config_dict)

    if modified_incar is not None:
        assert isinstance(modified_incar, dict), 'Invalid format for INCAR modification. Please provide a dictionnary with INCAR tags as keys and respective value as values'

    else:
        modified_incar = _MITRelaxSet_INCAR_corrections(structure.num_sites)

    if modified_kpoints is not None:
        assert isinstance(modified_kpoints, dict), 'Invalid format for KPOINTS modification. Please provide a dictionnary with KPOINTS PyMatGen supported mode as key and correct mesh definition as value'

    if modified_potcar is not None:
        assert isinstance(modified_potcar, dict), 'Invalid format for POTCAR modification. Please provide a dictionnary with element symbols as keys and a dictionnary containing hash and PP symbol as values'

    return MITRelaxSet(structure, user_incar_settings=modified_incar, user_kpoints_settings=modified_kpoints, user_potcar_settings=modified_potcar) 

def vasp_input_files_generator(structure: Structure, path: str, kwargs: Optional[Dict] = None):
    if kwargs is None:
        vasp_input_files_settings(structure).write_input(path)
    else:
        vasp_input_files_settings(structure, **kwargs).write_input(path)