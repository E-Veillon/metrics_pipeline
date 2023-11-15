#!/usr/bin/python

"""
A script using VASP DFT calculations to determine the relative stability of given structures in order to keep synthesizable structure and discard others.
"""

from typing import List
from itertools import count
from datetime import datetime
from argparse import ArgumentParser
from tqdm.contrib.concurrent import process_map
from pymatgen.core.structure import Structure


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
        default='./Vasp_input_sets/',
        help='Directory to write VASP input files in (created if it does not exist).',
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

    structures = read_cif(
        filename=args.datafile, 
        workers=args.workers, 
        keep_rare_gases=args.keep_rare_gases
        )
    nbr_struct      = len(structures)
    chunksize       = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
    generic_path    = args.path

    def feed_args(structures: List[Structure], path: str):
        if not path.endswith('/'):
            path += '/'
        struct_counter = count(0)
        return [(struct, path + f"{next(struct_counter)}/") for struct in structures]
    
    process_map(
        vasp_input_files_generator, 
        feed_args(structures, generic_path),  
        max_workers=args.workers, 
        chunksize=chunksize
    )

    end = datetime.now()
    print(end-start)

if __name__ == '__main__':
    main()
