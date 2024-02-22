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
from screening_pipeline.utils.typing import PMGStaticSet
from screening_pipeline.utils.vasp_io import vasp_batch_launch, vasp_launcher, batch_extract_vasp_data, \
                                             delta_sol_inputs_init
from screening_pipeline.utils.data_process import batch_calculate_delta_sol_band_gaps

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
        Missing steps to complete this script: None
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
        default='band_gap_passed.txt', 
        help='''Defines a file whose presence in a structure directory means it passed this
        screening step successfully and can be kept for further calculations.''',  
        metavar='accept_file.txt'
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
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to spawn for parallelized steps.',
        metavar='int',
    )
    parser.add_argument(
        '-t', 
        '--task_index', 
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
    workers        = args.workers
    struct_idx     = args.task_index


    # MAIN BLOCK

    # Extract relevant previous VASP outputs
    structs_data = batch_extract_vasp_data(
        method='delta_sol', 
        base_dir=input_dir, 
        ignore_file=ignore_file, 
        workers=workers
    )

    if struct_idx is None: # Use tqdm.contrib.concurrent.process_map() to do the calculation

        inputs_data = delta_sol_inputs_init(
            structs_data=structs_data, 
            preset=preset, 
            user_corrections=user_settings or None, 
            with_uncertainties=args.with_uncertainties
        )

        # Launch static calculations
        vasp_batch_launch(
            vasp_exe=exe_path, 
            inputs_data=inputs_data, 
            base_dir=outdir, 
            workers=workers
        )

        # Extract resulting energies
        bg_structs_data = batch_extract_vasp_data(
            method='delta_sol', 
            base_dir=outdir, 
            workers=workers
        )


    elif struct_idx >= 0: # Use the job array to parallelize the calculation

        try: struct_name = next(filter(
            lambda key: key.startswith(f'{struct_idx}_'), 
            structs_data.keys()
            ))
        except StopIteration:
            raise ValueError(f'Provided "task_index" arg is out of the range of indexed structures.')
        
        structs_data = {struct_name: structs_data.get(struct_name)}

        inputs_data = delta_sol_inputs_init(
            structs_data=structs_data, 
            preset=preset, 
            user_corrections=user_settings or None, 
            with_uncertainties=args.with_uncertainties
        )

        for name, vasp_input in inputs_data.items():
            run_path = os.path.join(outdir, name)
            vasp_launcher(vasp_exe=exe_path, path=run_path, vasp_input=vasp_input)

        bg_structs_data = batch_extract_vasp_data(
            method='delta_sol', 
            base_dir=outdir, 
            structs_names=list(inputs_data.keys()), 
            workers=workers
        )

    else: raise AssertionError('"task_index" arg must be positive or zero.')

    E_band_gaps = batch_calculate_delta_sol_band_gaps(
        structs_data, bg_structs_data, args.with_uncertainties, workers
    )

    good_bg_structs = list(filter(
        lambda tup: min(valid_interval) <= tup[1] <= max(valid_interval), 
        list(E_band_gaps.items())
    ))

    bad_bg_structs = list(filter(
        lambda tup: tup[1] < min(valid_interval) or tup[1] > max(valid_interval), 
        list(E_band_gaps.items())
    ))

    # Keep good structures
    for struct in good_bg_structs:
        name         = struct[0]
        E_band_gap   = round(struct[1], 6)
        name_plus    = '_'.join((name, 'best', 'plus'))
        name_minus   = '_'.join((name, 'best', 'minus'))
        accept_msg   = [
            "BAND GAP TEST PASSED", 
            f"Δ-Sol band gap was estimated to {E_band_gap} eV, which is inside the interval [{min(valid_interval)}, {max(valid_interval)}].", 
            "Therefore, it is suitable for wanted application, and should be considered for further screening steps."
        ]
        
        if args.with_uncertainties:
            E_band_gap_min = round(struct[2], 6)
            E_band_gap_max = round(struct[3], 6)
            accept_msg += [
                "\nUncertainty interval (does not affect acception or rejection):", 
                f"Band Gap minimum = {E_band_gap_min} eV", 
                f"Band Gap maximum = {E_band_gap_max} eV"
            ]

        accept_msg = "\n".join(accept_msg)
        input_path = os.path.join(input_dir, name, accept_file)
        path_plus  = os.path.join(outdir, name_plus, accept_file)
        path_minus = os.path.join(outdir, name_minus, accept_file)

        with open(input_path, 'wt') as f1, open(path_plus, 'wt') as f2, open(path_minus, 'wt') as f3:
            f1.write(accept_msg)
            f2.write(accept_msg)
            f3.write(accept_msg)

    # Reject unsuitable structures
    for struct in bad_bg_structs:
        name         = struct[0]
        E_band_gap   = round(struct[1], 6)
        name_plus    = '_'.join((name, 'best', 'plus'))
        name_minus   = '_'.join((name, 'best', 'minus'))
        reject_msg   = [
            "BAND GAP REJECTION", 
            f"Δ-Sol band gap was estimated to {E_band_gap} eV, which is not inside the interval [{min(valid_interval)}, {max(valid_interval)}].",
            "Therefore, it is not suitable for wanted application, and should not be considered in further screening steps."
        ]

        if args.with_uncertainties:
            E_band_gap_min = round(struct[2], 6)
            E_band_gap_max = round(struct[3], 6)
            reject_msg += [
                "\nUncertainty interval (does not affect acception or rejection):", 
                f"Band Gap minimum = {E_band_gap_min} eV", 
                f"Band Gap maximum = {E_band_gap_max} eV"
            ]

        reject_msg = "\n".join(reject_msg)
        input_path = os.path.join(input_dir, name, ignore_file)
        path_plus  = os.path.join(outdir, name_plus, ignore_file)
        path_minus = os.path.join(outdir, name_minus, ignore_file)

        with open(input_path, 'wt') as f1, open(path_plus, 'wt') as f2, open(path_minus, 'wt') as f3:
            f1.write(reject_msg)
            f2.write(reject_msg)
            f3.write(reject_msg)

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()
