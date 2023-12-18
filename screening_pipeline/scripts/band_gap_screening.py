#!/usr/bin/python
'''
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (reference 32 in screening_pipeline/Bibliography)
'''

########################################
# SYSTEM I/O MODULES

from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from pathlib import Path

########################################
# PYTHON MATERIALS GENOMICS PACKAGE


########################################
# LOCAL MODULES

# from screening_pipeline.utils.cif_io import read_cif
from screening_pipeline.utils.periodic_table import get_delta_sol_el_ratio
from screening_pipeline.utils.vasp_io import vasp_static_settings, vasp_batch_launch, \
                                             batch_extract_vasp_data, chgcar_density_switch, \
                                             delta_sol_inputs_init, batch_calculate_delta_sol_band_gaps, \
                                             filter_by_band_gap

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:
    
    allowed_presets = [
        'MPStaticSet', 
        'MatPESStaticSet', 
        'MPScanStaticSet'
    ]

    assert Path(args.input_dir).is_dir(), \
    f'{args.input_dir}: No directory found.'

    assert Path(args.output).is_dir(), \
    f'{args.output}: No directory found.'

    assert args.method in allowed_presets, \
    f'Provided static preset must be one of the following:\n \
    {allowed_presets}'

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
        Missing steps to complete this script:
            - Set an option to enable incertainty calculation on the band gap using N*min and N*max
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
        '-m', 
        '--method', 
        type=str, 
        default='MPStaticSet', 
        help='''The pymatgen preset to use for VASP static calculations.
                More info on possible presets in pymatgen documentation:
                https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.''', 
        metavar='StaticSet'
    )
    parser.add_argument(
        '-f', 
        '--functional', 
        type=str, 
        default='PBE', 
        help='DFT functional to use for Δ-Sol N* parameter initialization.', 
        metavar='str'
    )
    parser.add_argument(
        '-n', 
        '--n_star_type', 
        type=str, 
        default='BEST', 
        help='Type of N* parameter to initialize for Δ-Sol method.', 
        metavar='str'
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
    args: Namespace = parser.parse_args()

    assert_args(args)
    
    input_dir   = Path(args.input_dir)
    outdir      = Path(args.output)
    preset      = args.method
    functional  = args.functional
    n_star_type = args.n_star_type
    ignore_file = args.ignore
    workers     = args.workers


    # MAIN BLOCK

    # Extract relevant previous VASP outputs
    structs_data = batch_extract_vasp_data(
        method='delta_sol', 
        base_dir=input_dir, 
        ignore_file=ignore_file, 
        workers=workers
    )

    inputs_list, subdirs_list = delta_sol_inputs_init(
        structs_data=structs_data, 
        preset=preset, 
        dft_functional=functional, 
        n_star_type=n_star_type
    )

    # Launch static calculations
    vasp_batch_launch(
        vasp_inputs=inputs_list, 
        base_dir=outdir, 
        subdir_names=subdirs_list, 
        worker=workers
    )

    # Extract resulting energies
    bg_structs_data = batch_extract_vasp_data(
        method='delta_sol', 
        base_dir=outdir, 
        workers=workers
    )

    E_band_gaps = batch_calculate_delta_sol_band_gaps(
        structs_data, bg_structs_data, workers
    )

    # Reject unsuitable structures
    filter_by_band_gap(
        E_band_gaps, 
        valid_interval=[1.3, 3.6], 
        base_dir=outdir, 
        ignore_file=ignore_file
    )

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()