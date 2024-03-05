#!/usr/bin/python
'''
A script that writes and performs VASP relaxation on structure data read from CIF files,
using pymatgen as a setting interface between raw data and VASP.
The relaxation results then may be used in other scripts for material properties analysis.
'''

########################################
# SYSTEM I/O MODULES

import os
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter

########################################
# OPTIMIZATION MODULES

from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import SiteCollection

########################################
# LOCAL MODULES

from screening_pipeline.utils.utils import _yaml_loader
from screening_pipeline.utils.typing import PMGRelaxSet
from screening_pipeline.utils.cif_io import read_cif
from screening_pipeline.utils.vasp_io import vasp_relaxation_settings, vasp_batch_launch, vasp_launcher

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:
    
    assert os.path.exists(args.filename), (
    f'{args.filename}: path to input file not found.'
    )
    assert os.path.isfile(args.filename), (
    f'{args.filename} found but it is not a file.'
    )
    assert args.filename.endswith('.cif'), (
    'Input structure data must be in CIF format.'
    )
    assert args.executable_path.startswith("vasp") or os.path.exists(args.executable_path), (
    f'{args.executable_path}: executable file not found.'
    )
    assert os.path.exists(args.output), (
    f'{args.output}: path to output directory not found.'
    )
    assert os.path.isdir(args.output), (
    f'{args.output} found but it is not a directory.'
    )
    assert args.preset in PMGRelaxSet, (
    'Provided relaxation preset must be one of the following:\n'
    f'{PMGRelaxSet}'
    )
    assert os.path.isfile(args.user_settings), (
    f'{args.user_settings}: file not found.'
    )
    assert args.user_settings.endswith('.yaml'), (
    'user settings file must be of .yaml format.'
    )
    assert args.workers >= 1, (
    '"workers" arg must be strictly positive.'
    )
    if args.task_index is not None:
        assert args.task_index >= 0 or args.task_index is None, (
        "'task_index' arg must be positive or zero."
        )

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
        metavar='input_file.cif'
    )
    parser.add_argument(
        'executable_path', 
        type=str,
        help='Path to the VASP executable.', 
        metavar='/path/to/vasp'
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
        '-p', 
        '--preset', 
        type=str, 
        default='MPRelaxSet', 
        help='''The pymatgen preset to use for VASP relaxation.
                More info on possible presets in pymatgen documentation:
                https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.''', 
        metavar='RelaxSet'
    )
    parser.add_argument(
        '-u', 
        '--user-settings', 
        type=str, 
        default='user_settings.yaml', 
        help='Path to the .yaml file containing tags overrides to put over the PMG preset.', 
        metavar='file.yaml', 
        dest='user_settings'
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to spawn.',
        metavar='int',
    )
    parser.add_argument(
        '-t', 
        '--task_index', 
        type=int, 
        help='If a job array is used, provide here the structure index to treat according to task IDs\n \
            (e.g. if task ID 0 treats structure 0 and so on, just provide the task ID).'
    )

    args: Namespace = parser.parse_args()

    assert_args(args)
    
    input_file    = args.filename
    exe_path      = args.executable_path
    outdir        = args.output
    preset        = args.preset
    user_settings = _yaml_loader(args.user_settings)
    workers       = args.workers
    struct_idx    = args.task_index


    # MAIN BLOCK

    # Convert CIF data into Structure objects
    structures, *_ = read_cif(
        filename=input_file, 
        keep_rare_gases=True, # Avoid calling rare gaz screening function
        keep_rare_earths=True # Avoid calling rare earth screening function
    )

    if struct_idx is None: # Use tqdm.contrib.concurrent.process_map() to do the calculation
        
        #TODO: Revoir cette partie puisqu'elle ne fonctionne que pour une structure à la fois

        # Setup parallel processing
        nbr_struct      = len(structures)
        chunksize       = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
        vasp_setup      = partial(
            vasp_relaxation_settings, 
            preset=preset, 
            user_corrections=user_settings
        )
        dir_names_list  = [f'{idx}_{structure.composition.reduced_formula}' for idx, structure in enumerate(structures)]

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
            vasp_exe=exe_path, 
            vasp_inputs=vasp_inputs, 
            base_dir=outdir, 
            subdir_names=dir_names_list, 
            workers=workers
        )
    
    elif struct_idx >= 0: # Use the job array to parallelize the calculation

        structure: SiteCollection = structures[struct_idx]
        dir_name                  = f'{struct_idx}_{structure.composition.reduced_formula}'
        
        vasp_input = vasp_relaxation_settings(
            structure=structure, 
            preset=preset, 
            user_corrections=user_settings
        )

        vasp_launcher(vasp_exe=exe_path, vasp_input=vasp_input, path=os.path.join(outdir, dir_name))
    
    else: raise AssertionError('"task_index" arg must be positive or zero.')

    stop = datetime.now()
    print(f'Elapsed time: {stop-start}')


if __name__=='__main__':
    main()
