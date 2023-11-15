#!/usr/bin/python

"""
A script using VASP DFT calculations to determine the relative stability of given structures in order to keep synthesizable structure and discard others.
"""

from datetime import datetime
from argparse import ArgumentParser


def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    name     = 'stability_screening.py'
    desc     = 'A script using VASP DFT calculations to determine the relative stability of given structures in order to keep synthesizable structure and discard others.'
    footnote = 'A first try to writing cleaner code (unfinished).'

    parser = ArgumentParser(prog=name, description=desc, epilog=footnote)
  
    parser.add_argument(
        'datafile',
        type=str,
        help='The file conaining structure data to write VASP input files for. CIF format only.',
        metavar='str'
    )
    parser.add_argument(
        '-p', '--path',
        type=str,
        default='.',
        help='Directory to write VASP input files in.',
        metavar='str'
    )
    parser.add_argument(
        '-w', '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create',
        metavar='int'
    )
    parser.add_argument(
        "--keep-rare-gases",
        action="store_true",
        help="Pass this flag to disable automatic elimination of structures containing rare gases"
    )

    args = parser.parse_args()

    # MAIN BLOCK

    from screening_pipeline.utils import read_cif, vasp_input_files_generator

    structs         = read_cif(filename=args.datafile, workers=args.workers, keep_rare_gases=args.keep_rare_gases)
    path            = args.path
    Natoms          = struct.num_sites #number of atoms in the unit cell (depending on the structure)
    corrected_EDIFF = float(5e-5)*Natoms
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
            'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4.0, 'Fe': 4.0, #'Cu': 4 -> 4.0
            'Mn': 3.9, 'Mo': 3.5, 'Nb': 1.5, 'Ni': 6.0, 'Re': 2.0, #'Mo': 4.38 -> 3.5, 'Ni': 6 -> 6.0, 'Re': 2 -> 2.0
            'Ta': 2.0, 'V': 3.1, 'W': 4.0                          #'Ta': 2 -> 2.0
        }, 
        'O': {
            'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4.0, 'Fe': 4.0, #'Cu': 4 -> 4.0
            'Mn': 3.9, 'Mo': 3.5, 'Nb': 1.5, 'Ni': 6.0, 'Re': 2.0, #'Mo': 4.38 -> 3.5, 'Ni': 6 -> 6.0, 'Re': 2 -> 2.0
            'Ta': 2.0, 'V': 3.1, 'W': 4.0                          #'Ta': 2 -> 2.0
        }, 
        'S': {
            'Fe': 1.9, 'Mn': 2.5
        }}
    corrected_INCAR = {
        "EDIFF": corrected_EDIFF,
        "ENCUT": corrected_ENCUT,
        "LDAUL": corrected_LDAUL,
        "LDAUU": corrected_LDAUU, 
        "LMAXMIX": 4 #Necessary to get reliable results with GGA + U framework
    }
    Input_dict = vasp_input_files_generator(struct, modified_incar=corrected_INCAR)
    Input_dict.write_input(path)
    end = datetime.now()
    print(end-start)

if __name__ == '__main__':
    main()
