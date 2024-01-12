#!/usr/bin/python
'''
A script that writes and performs VASP relaxation on structure data read from CIF files,
using pymatgen as a setting interface between raw data and VASP.
The relaxation results then may be used in other scripts for material properties analysis.
'''

########################################
# SYSTEM I/O MODULES

from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from pathlib import Path

########################################
# OPTIMIZATION MODULES

from itertools import count
from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# LOCAL MODULES

from screening_pipeline.utils.cif_io import read_cif
from screening_pipeline.utils.vasp_io import vasp_relaxation_settings, vasp_batch_launch

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:
    
    assert Path(args.filename).exists(), \
    f'{args.filename}: path to input file not found.'
    
    assert Path(args.filename).is_file(), \
    f'{args.filename} found but it is not a file.'

    assert args.filename.endswith('.cif'), \
    'Input structure data must be in CIF format.'

    assert Path(args.output).exists(), \
    f'{args.output}: path to output directory not found.'

    assert Path(args.output).is_dir(), \
    f'{args.output} found but it is not a directory.'

    allowed_presets = {
        'MITRelaxSet', 
        'MPRelaxSet', 
        'MPScanRelaxSet', 
        'MPHSERelaxSet', 
        'MPMetalRelaxSet', 
        'MVLRelax52Set', 
        'MVLScanRelaxSet'
    }

    assert args.method in allowed_presets, \
    f'Provided relaxation preset must be one of the following:\n \
    {allowed_presets}'

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = 'vasp_relaxation.py'
    prog_description = '''
    A script that writes and performs VASP relaxation on structure data read from CIF files,
    using pymatgen as a setting interface between raw data and VASP.
    The relaxation results then may be used in other scripts for material properties analysis.
    '''
    prog_missing_steps = '''
        Missing steps to complete this script: None
        '''
    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_description, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        'filename',
        type=str,
        help='Path to the CIF file containing structure data to read.', 
        metavar='input_file'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default='./',
        help='''Path to the output directory where VASP files will be written.
                A subdirectory will be created in output directory
                for each structure found in file.cif.''', 
        metavar='outdir'
    )
    parser.add_argument(
        '-m', 
        '--method', 
        type=str, 
        default='MITRelaxSet', 
        help='''The pymatgen preset to use for VASP relaxation.
                More info on possible presets in pymatgen documentation:
                https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.''', 
        metavar='RelaxSet'
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create',
        metavar='int',
    )
    args: Namespace = parser.parse_args()

    assert_args(args)
    
    input_file = Path(args.filename)
    outdir     = Path(args.output)
    preset     = args.method
    workers    = args.workers


    # MAIN BLOCK

    # Convert CIF data into Structure objects
    structures, *_ = read_cif(filename=input_file)
    
    # Setup parallel processing
    nbr_struct      = len(structures)
    chunksize       = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
    vasp_setup      = partial(vasp_relaxation_settings, preset=preset)
    dir_names_list  = []
    
    for idx, structure in enumerate(structures):
        struct_dir_name = f'{idx}_{structure.formula}'
        dir_names_list.append(struct_dir_name)

    # Write VaspInput objects from structures and chosen preset
    vasp_inputs = list(process_map(
        vasp_setup, 
        structures, 
        max_workers=workers, 
        chunksize=chunksize, 
        desc='Writing VASP input files'
    ))

    # Use written VaspInput objects to write input files and run VASP
    vasp_batch_launch(
        vasp_input=vasp_inputs, 
        base_dir=outdir, 
        subdir_names=dir_names_list, 
        workers=workers
    )

    stop = datetime.now()
    print(f'Elapsed time: {stop-start}')


if __name__=='__main__':
    main()