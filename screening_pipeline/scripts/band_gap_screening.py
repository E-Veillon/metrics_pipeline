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

#from pymatgen.io.vasp.outputs import Oszicar, Chgcar

########################################
# LOCAL MODULES

# from screening_pipeline.utils.cif_io import read_cif
from screening_pipeline.utils.periodic_table import get_delta_sol_el_ratio
from screening_pipeline.utils.vasp_io import vasp_static_settings, vasp_batch_launch, \
                                             batch_extract_vasp_data, chgcar_density_switch

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
            - Add necessary in-line command arguments
            - Calculate N0, number of valence electrons in unit cell - OK
            - Choose right N*best (PBE[sp] or PBE[spd]) - OK
            - n = N0/N*best - OK
            - Calculate E(N0), E(N0 + n), E(N0 - n) with VASP (static ionically)
            - Calculate E(gap) = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n
            - For photocatalysts article: keep materials with 1.3 < E(gap) < 3.6 eV by default
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
    ignore_file = args.ignore
    workers     = args.workers


    # MAIN BLOCK

    # Extract relevant previous VASP outputs
    structs_data = batch_extract_vasp_data(
        base_dir=input_dir, 
        ignore_file=ignore_file, 
        workers=workers
    )

    # TODO: Le bloc ci-dessous reste à paralléliser
    inputs_list  = []
    subdirs_list = []
    for name, data in structs_data.items():
        # Bloc de data-splitting
        structure = data['structure']
        chgcar    = data['CHGCAR']
        delta     = get_delta_sol_el_ratio(structure)
        data['n_ratio'] = delta
        data['CHGCAR_plus'], data['CHGCAR_minus'] = chgcar_density_switch(chgcar, delta)
        # Bloc de préparation des calculs
        run_plus  = vasp_static_settings(structure)
        run_plus.update({'CHGCAR': data['CHGCAR_plus']})
        run_plus_path  = Path('_'.join(name , 'plus'))
        run_minus = vasp_static_settings(structure)
        run_minus.update({'CHGCAR': data['CHGCAR_minus']})
        run_minus_path = Path('_'.join(name , 'minus'))
        inputs_list.extend([run_plus, run_minus])
        subdirs_list.extend([run_plus_path, run_minus_path])
    
    # Launch static calculations
    vasp_batch_launch(
        vasp_inputs=inputs_list, 
        base_dir=outdir, 
        subdir_names=subdirs_list, 
        worker=workers
    )

    # Extract resulting energies
    bg_structs_data = batch_extract_vasp_data(
        base_dir=outdir, 
        workers=workers
    )

    final_energies = {name: data['final_energy'] for name, data in bg_structs_data.items()}

    # E_FG = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n -> Δ-Sol band gap 
    # (Ref 32 in screening_pipeline/Bibliography))
    # TODO: Reste à paralléliser
    for name, data in structs_data.items():
        name_plus    = '_'.join(name, 'plus')
        name_minus   = '_'.join(name, 'minus')
        E_N0         = data['final_energy']
        E_N0_plus_n  = final_energies[name_plus]['final_energy']
        E_N0_minus_n = final_energies[name_minus]['final_energy']
        n            = data['n_ratio']
        E_band_gap   = (E_N0_plus_n + E_N0_minus_n - 2*E_N0)/n
        bg_too_small = E_band_gap < 1.3
        bg_too_big   = E_band_gap > 3.6

        # Reject unsuitable structures
        if bg_too_small or bg_too_big:
            reject_str = f'Δ-Sol band gap was estimated to {E_band_gap} eV, which is not inside the interval [1.3, 3.6].\n \
            Therefore, it is not suitable for wanted application and should not be considered in further screening steps.'
            
            path_plus  = Path(outdir / name_plus / ignore_file)
            path_plus.touch()
            path_plus.write_text(reject_str)
            
            path_minus = Path(outdir / name_minus / ignore_file)
            path_minus.touch()
            path_minus.write_text(reject_str)

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()