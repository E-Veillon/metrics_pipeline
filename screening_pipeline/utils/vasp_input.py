#!/usr/bin/python

from typing import Optional
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import DictSet, MITRelaxSet 

def vasp_input_files_generator(structure: Structure, /, *, use_mit_set: bool = True, config_dict: Optional[dict] = None, modified_incar: Optional[dict] = None, modified_kpoints: Optional[dict] = None, modified_potcar: Optional[dict] = None):
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

            modified_incar (dict):      User INCAR settings. It allows to override some of the standards INCAR tags if necessary.
                                        Defaults to None.

            modifieed_kpoints (dict):   User KPOINTS settings. It allows to override the standard Kpoints mesh if necessary..
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

    if modified_kpoints is not None:
        assert isinstance(modified_kpoints, dict), 'Invalid format for KPOINTS modification. Please provide a dictionnary with KPOINTS PyMatGen supported mode as key and correct mesh definition as value'

    if modified_potcar is not None:
        assert isinstance(modified_potcar, dict), 'Invalid format for POTCAR modification. Please provide a dictionnary with element symbols as keys and a dictionnary containing hash and PP symbol as values'

    return MITRelaxSet(structure, user_incar_settings=modified_incar, user_kpoints_settings=modified_kpoints, user_potcar_settings=modified_potcar) 
