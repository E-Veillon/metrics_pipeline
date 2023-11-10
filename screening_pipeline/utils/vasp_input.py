#!/usr/bin/python

from typing import Optional
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MITRelaxSet

def vasp_input_files_generator(structure: Structure, modified_incar: Optional[dict] = None, modified_kpoints: Optional[dict] = None, modified_potcar: Optional[dict] = None):
    '''
    Build a VASP input files generator with a standard preset corresponding to the MIT high throughput material screening project.
        
        Reference:
            A. Jain, G. Hautier, C.J. Moore, S.P. Ong,
            C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder,
            Computational Materials Science, 50, 2295-2310 (2011)

        Parameters:
            structure (Structure): The structure to write VASP input files for.

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

    if modified_incar is not None:
        assert isinstance(modified_incar, dict), 'Invalid format for INCAR modification. Please provide a dictionnary with INCAR tags as keys and respective value as values'

    if modified_kpoints is not None:
        assert isinstance(modified_kpoints, dict), 'Invalid format for KPOINTS modification. Please provide a dictionnary with KPOINTS PyMatGen supported mode as key and correct mesh definition as value'

    if modified_potcar is not None:
        assert isinstance(modified_potcar, dict), 'Invalid format for POTCAR modification. Please provide a dictionnary with element symbols as keys and a dictionnary containing hash and PP symbol as values'

    return MITRelaxSet(structure=structure, user_incar_settings=modified_incar, user_kpoints_settings=modified_kpoints, user_potcar_settings=modified_potcar) 
