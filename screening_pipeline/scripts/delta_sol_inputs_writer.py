#!/usr/bin/python
'''
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (reference 32 in screening_pipeline/Bibliography)
'''

########################################
# SYSTEM I/O MODULES

import os
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter

########################################
# PYTHON MATERIALS GENOMICS PACKAGE


########################################
# LOCAL MODULES

from screening_pipeline.utils.utils import _yaml_loader
from screening_pipeline.utils.custom_types import PMGStaticSet
from screening_pipeline.utils.vasp_io import (
    vasp_launcher, batch_extract_vasp_data, 
    delta_sol_inputs_init, delta_sol_calculation_init
)

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:

    assert os.path.isdir(args.input_dir), \
    f'{args.input_dir}: No directory found.'

    assert os.path.exists(args.executable_path), \
    f'{args.executable_path}: No such file found.'

    assert os.path.isdir(args.output), \
    f'{args.output}: No directory found.'

    assert args.preset in PMGStaticSet, \
    f'Provided static preset must be one of the following:\n \
    {PMGStaticSet}'

    assert os.path.isfile(args.user_settings), \
    f'{args.user_settings}: file not found.'

    assert args.user_settings.endswith('.yaml'), \
    'user settings file must be of .yaml format.'

    assert args.mini_maxi[0] >= 0.0 and args.mini_maxi[1] >= 0.0, \
    f"Acceptable band gap values must be positive or zero."

    assert args.mini_maxi[0] != args.mini_maxi[1], \
    f"Acceptable band gap values cannot have the same value."

    assert args.accept.endswith('.txt'), \
    'Structure acceptance file must be a plain text file type (.txt).'

    assert args.ignore.endswith('.txt'), \
    'Structure ignoring file must be a plain text file type (.txt).'

    assert args.workers >= 1, \
    'The number of workers cannot be negative or zero.'

    assert args.task_index >= 0 or args.task_index is None, \
    '"task_index" arg must be positive or zero.'

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK
    
    prog_name = 'band_gap_screening'
    prog_desc = '''
        A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

        Reference for Δ-Sol method:
            M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
            (reference 32 in screening_pipeline/Bibliography)
        '''
    prog_missing_steps = '''
        Missing steps to complete this script:
            - Set an option to enable incertainty calculation on the band gap using N*min and N*max - OK
            - This option should let user choose if they want to keep materials according to the uncertainty case:
                * If E(gap) is inside the goal but uncertainty gets out ?
                * If E(gap) is outside the goal but uncertainty gets in ?
                Possible options: keep, keep_aside, discard, auto (keep_aside if at least half the interval is in, else discard)
        '''
    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_desc, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        'input_dir',
        type=str,
        help='Base directory containing structure directories.', 
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
                for each structure processed.''', 
        metavar='outdir'
    )
    parser.add_argument(
        '-p', 
        '--preset', 
        type=str, 
        default='MPStaticSet', 
        help='''The pymatgen preset to use for VASP static calculations.
                More info on possible presets in pymatgen documentation:
                https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.''', 
        metavar='StaticSet'
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
        '-m',
        '--mini-maxi',
        nargs=2,
        type=float,
        default=[1.3, 3.6],
        help='''Acceptable interval of band gap values in eV (Defaults: %(default)s eV).''', 
        metavar='float float'
    )
    parser.add_argument(
        '-a', 
        '--accept', 
        type=str, 
        default=None, 
        help='''Defines a file whose presence in a structure directory means it passed this
        screening step successfully and can be kept for further calculations.''',  
        metavar='accept_file.txt'
    )
    parser.add_argument(
        '-i', 
        '--ignore', 
        type=str, 
        default=None, 
        help='''Defines a file whose presence in a structure directory means it did not pass
        previous screening steps and should not be used in this calculation. 
        This file will also be written in structure directories that did not pass this step.
        WARNING: 
        If the name of this file is overwritten, care must be taken that it is the same file
        throughout every used screening steps to make sure rejected structures don't go further.''', 
        metavar='ignore_file.txt'
    )
    parser.add_argument(
        "-s",
        "--summary",
        default="summary.json",
        help=(
            "Output file indicating calculation results for this step (json format).\n"
            "Also used by further steps to filter out structures that were rejected in previous steps."
        ),
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to spawn for parallelized steps.',
        metavar='int',
    )
    parser.add_argument(
        '-t', 
        '--task-index', 
        type=int, 
        help='If a job array is used, provide here the structure index to treat according to task IDs\n \
            (e.g. if task ID 0 treats structure 0 and so on, just provide the task ID).'
    )
    parser.add_argument(
        '--with-uncertainties', 
        action='store_true', 
        help='Pass this flag to enable computation of minimal and maximal Δ-Sol band gaps.\n\
              This will need one full VASP static total energy computation for each limit.'
    )

    args: Namespace = parser.parse_args()

    assert_args(args)

    input_dir      = args.input_dir
    exe_path       = args.executable_path
    outdir         = args.output
    preset         = args.preset
    user_settings  = _yaml_loader(args.user_settings, on_error='raise')
    valid_interval = sorted(args.mini_maxi)
    accept_file    = args.accept
    ignore_file    = args.ignore
    summary_file   = args.summary
    workers        = args.workers


    # MAIN BLOCK

    # Extract relevant previous VASP outputs
    structs_data = batch_extract_vasp_data(
        method='delta_sol', 
        base_dir=input_dir, 
        ignore_file=ignore_file, 
        path_to_summary=summary_file, 
        workers=workers
    )

    if 'task_index' in args: # Initialize input and run VASP on it
        task_index       = args.task_index
        tasks_per_struct = 7 if args.with_uncertainties else 3
        struct_idx       = task_index // tasks_per_struct
        calc_idx         = task_index % tasks_per_struct

        try:
            struct_name = next(filter(
                lambda key: key.startswith(f'{struct_idx}_'), 
                structs_data.keys()
            ))
        except StopIteration:
            raise ValueError(
                "Provided 'task_index' arg is out of the range of indexed structures, "
                "or the corresponding structure is already rejected."
            )
    
        input_data = delta_sol_calculation_init(
            structure=structs_data[struct_name].get('structure'), 
            calc_index=calc_idx, 
            preset=preset, 
            user_corrections=user_settings
        )
        vasp_launcher(vasp_exe=exe_path, path=outdir, vasp_input=input_data)

    elif not "task_index" in args: # Write inputs only
        structs_data = {struct_name: structs_data.get(struct_name)}

        inputs_data = delta_sol_inputs_init(
            structs_data=structs_data, 
            preset=preset, 
            user_corrections=user_settings or None, 
            with_uncertainties=args.with_uncertainties
        )

        for name, vasp_input in inputs_data.items():
            run_dir = os.path.join(outdir, name)
            vasp_input.write_input(output_dir=run_dir)

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()
