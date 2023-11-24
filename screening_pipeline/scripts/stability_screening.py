#!/usr/bin/python

'''
A script using VASP DFT calculations to determine the relative stability of given structures.
'''

from typing import List, Tuple
from itertools import count
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter
from tqdm.contrib.concurrent import process_map
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import VaspInput

def assert_args(args: Namespace):
    '''
    Input arguments verification.

    Args:
        args (NamedTuple): namespace of the parsed arguments.
    '''

    assert args.filename.endswith('.cif'), 'Input file must be in CIF format'
    assert args.workers >= 1, 'the number of workers cannot be negative or zero'

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name          = 'stability_screening'
    prog_description   = 'A script using VASP DFT calculations to determine the relative stability of given structures.'
    prog_missing_steps = '''
        Missing steps to complete the script:
            - Determine the critical formation reaction of the structure
            - (Relax reference structures for calculations consistency)
            - Compare structure energy with the sum of reference structure energies
            - Aknowledge what to do with unstable structures
        '''
    helper_format = RawDescriptionHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_description, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )
  
    parser.add_argument(
        'filename',
        type=str,
        help='The file conaining structure data to calculate relative instability from. CIF format only.'
    )
    parser.add_argument(
        '-p', '--path',
        type=str,
        default='./Vasp_input_sets/',
        help='Directory to write VASP input files in (created if it does not exist) (default = ./Vasp_input_sets/).',
        metavar='str'
    )
    parser.add_argument(
        '-w', '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create (default = 1).',
        metavar='int'
    )
    parser.add_argument(
        '--keep-rare-gases',
        action='store_true',
        help='Pass this flag to disable automatic elimination of structures containing rare gases.'
    )

    args: Namespace = parser.parse_args()

    assert_args(args)


    # MAIN BLOCK

    from screening_pipeline.utils import read_cif, vasp_input_files_settings, vasp_launcher

    structures = read_cif(
        filename=args.filename, 
        workers=args.workers, 
        keep_rare_gases=args.keep_rare_gases
    )

    nbr_struct: int   = len(structures)
    chunksize: int    = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
    generic_path: str = args.path

    def feed_args(vasp_inputs: List[VaspInput], path: str) -> List[Tuple]:
        if not path.endswith('/'):
            path += '/'
        calc_counter = count(0)
        return [(input, path + f'{next(calc_counter)}/') for input in vasp_inputs]
    
    vasp_input_sets: List[VaspInput] = process_map(
        vasp_input_files_settings, 
        structures, 
        max_workers=args.workers, 
        chunksize=chunksize
    )

    process_map(
        vasp_launcher, 
        feed_args(vasp_input_sets, generic_path), 
        max_workers=args.workers, 
        chunksize=chunksize
    )

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')

if __name__ == '__main__':
    main()
