#!/usr/bin/python
'''
A script using previous VASP relaxations to compute the relative stability of given structures.

Algorithm:

    1/  Group structures by dimension and composition. - OK
        
    2/  Get reference structures from a dataset and add them to the pool
        with a specific attribute. - Not done yet

    3/  Build one convex hull per structure group 
        (initialize ref elements from group composition) - OK
        
    4/  Inside each convex hull, compute ΔH for each generated entry. - OK

    5/  Reject structures that have a ΔH above given threshold value. - OK

Possible alternatives:

    Comparing to known references:
        1 - For each group, search in the dataset all corresponding structures.
        2 - Build the phase diagram corresponding to dataset structures.
        3 - Compute ΔH for all generated structures of the group.
        4 - Reject structures too much above reference convex hull.

    Opti:   Store some common phase spaces in a local util file.

    Opti:   Match and ignore structures equivalents to references with StructureMatcher.
'''

########################################
# SYSTEM I/O MODULES

from typing import List, Tuple
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from pathlib import Path

########################################
# OPTIMIZATION MODULES

from itertools import count
from tqdm.contrib.concurrent import process_map

########################################
# PYTON MATERIALS GENOMICS PACKAGE

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import VaspInput

########################################
# LOCAL MODULES

from screening_pipeline.utils.vasp_io import batch_extract_vasp_data
from screening_pipeline.utils.data_process import batch_calculate_instability_energies

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:

    assert Path(args.input_dir).is_dir(), \
    f'{args.input_dir}: No directory found.'

    assert args.ignore.endswith('.txt'), \
    'Structure ignoring file must be a plain text file type (.txt).'

    assert args.eliminate >= 1e-8, \
    '''Instability energy elimination criterion must be strictly positive.
       Moreover, any value below 1.10-5 meV/atom is too low and not supported.'''

    assert args.workers >= 1, \
    'The number of workers cannot be negative or zero.'

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name          = 'stability_screening'
    prog_description   = 'A script using previous VASP relaxations to compute the relative \
                          stability of given structures.'
    prog_missing_steps = '''
        Missing steps to complete the script:
            - Determine the critical formation reaction of the structure
            - (Relax reference structures for calculations consistency)
            - Compare structure energy with the sum of reference structure energies
            - Aknowledge what to do with unstable structures
        '''
    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_description, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )
  
    parser.add_argument(
        'input_dir',
        type=str,
        help='Base directory containing structure directories.', 
    )
    parser.add_argument(
        '-i', 
        '--ignore', 
        type=str, 
        default='rejected.txt', 
        help='''Defines a file whose presence in a structure directory means it did not pass
        previous screening steps and should not be used in this calculation. 
        This file will also be written in structure directories that did not pass this step.
        WARNING: 
        If the name of this file is overwritten, care must be taken that it is the same file
        throughout every used screening steps to make sure rejected structures don't go further.''', 
        metavar='ignore_file.txt'
    )
    parser.add_argument(
        '-e',
        '--eliminate',
        type=float,
        default=0.036,
        help=f'''Maximum value of ΔH (in eV/atom) above which structures are considered too unstable and rejected.
                 Defaults to 36 meV/atom, as used in the following paper, and seems fairly strict: 
                 Y. Wu, P. Lazic, G. Hautier, K. Persson, and G. Ceder, 
                 First principles high throughput screening of oxynitrides for water-splitting photocatalysts, 
                 Energy & Environmental Science 6, no. 1 (2012) 157''',
        metavar='float',
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to spawn for parallelized steps.',
        metavar='int',
    )
    args: Namespace = parser.parse_args()

    assert_args(args)
    
    input_dir     = Path(args.input_dir)
    ignore_file   = args.ignore
    delta_H_limit = round(args.eliminate, 8)
    workers       = args.workers


    # MAIN BLOCK

    structs_data = batch_extract_vasp_data(
        method='convex_hull', 
        base_dir=input_dir, 
        ignore_file=ignore_file, 
        workers=workers
    )

    structs_data = batch_calculate_instability_energies(structs_data, workers=workers)

    unstable_structs = list(filter(
        lambda _, data: data['delta_H'] > delta_H_limit, 
        structs_data.items()
    ))

    for struct in unstable_structs:
        name = struct[0]
        data = struct[1]
        reject_msg = f'''
                    Instability energy for this structure is estimated at {data['delta_H']} eV/atom, 
                    which is above the fixed instability limit of {delta_H_limit} eV/atom.
                    Therefore, it is considered not suitable for wanted application, 
                    and should not be considered in further screening steps.
                    '''
        reject_file_path = Path(input_dir / name / ignore_file)
        reject_file_path.touch()
        reject_file_path.write_text(reject_msg)

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')

if __name__ == '__main__':
    main()
