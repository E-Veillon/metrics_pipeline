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
    #ENMAX           = "Max cutoff energy automatically set in the POTCAR file" #Vérifier ENMAX une fois POTCAR généré une première fois pour entrer une valeur directe pour ENCUT (on veut ENCUT = 1.3*ENMAX)
    corrected_incar = dict("EDIFF": (5e-5)*Natoms, "LMAXMIX": 4, "ENCUT": 520, "LDAUL": {"S": {"Mn": 2}}, "LDAUU": {"F": {"Mo": 3.5}, "O": {"Mo": 3.5}})
    Input_dict = vasp_input_files_generator(struct, corrected_incar)
    Input_dict.write_input(path)
    end = datetime.now()
    print(end-start)

if __name__ == '__main__':
    main()
