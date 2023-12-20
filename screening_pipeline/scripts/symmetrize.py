#!/usr/bin/python
'''
A script that calculates symmetry spacegroup for structures in a CIF file using pymatgen.
'''

########################################
# SYSTEM I/O MODULES

from datetime import datetime
from argparse import ArgumentParser, RawTextHelpFormatter, Namespace

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace):
    '''
    Input arguments verification.

    Parameters:
        args (Namespace): namespace of the parsed arguments.
    '''

    print(' - INPUT ARGUMENTS - ')
    print(f'filename: {args.filename}')
    print(f'output: {args.output}')

    assert (
        args.filename.endswith('.cif')
        and args.output.endswith('.cif')
    ), 'some arguments formats are not supported, please only use CIF format'

    print(f'Tolerance in relative atomic distances check: {args.valid_tol}')

    assert args.valid_tol >= 0.0, \
    'relative atomic positions tolerance cannot be negative.'

    print(f'Fractional coordinates precision: {args.precision}')

    assert (
        0.0 <= args.precision <= 0.5
    ), 'Fractional coordinates precision must be between 0.0 and 0.5 to retain some reliability'

    print(f'Angles tolerance: {args.angleprec} degrees')

    assert (
        0.0 <= args.angleprec <= 20.0
    ), 'Angles tolerance must be between 0.0 and 20.0 degrees to retain some reliability'

    print(f'workers: {args.workers}')

    assert args.workers >= 1, 'the number of workers cannot be negative or zero'

    print(' ')
    print(' - OPTIONAL FLAGS - ')
    print(f'keep_rare_gases: {args.keep_rare_gases}')
    print(f'keep_equivalent: {args.keep_equivalent}')
    print(' ')
    print('------------------------------')
    print(' ')

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = 'symmetrize'
    prog_description = 'A script that calculates symmetry spacegroup for structures in a CIF file using pymatgen.'
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
        help='name of the file containing the input structures.'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default='[filename]_out.cif',
        help='name of output file containing filtered and symmetrized structures.',
    )
    parser.add_argument(
        '-d', 
        '--distance_tol', 
        type=float, 
        default=0.0, 
        help='''Tolerance for checking atoms relative positions in Angstroms.
                Structures containing atoms that are closer than this value will be discarded.
                The default value of 0.0 disables this feature.''', 
        metavar='float', 
        dest='valid_tol'
    )
    parser.add_argument(
        '-p',
        '--precision',
        type=float,
        default=0.01,
        help='Fractional coordinates tolerance for symmetry finding.',
        metavar='float',
    )
    parser.add_argument(
        '-a',
        '--angleprec',
        type=float,
        default=5.0,
        help='Angle tolerance for symmetry finding in degrees.',
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
    parser.add_argument(
        '--keep_rare_gases',
        action='store_true',
        help='A flag to disable automatic elimination of structures containing rare gases.'
    )
    parser.add_argument(
        '--keep_rare_earths',
        action='store_true',
        help='A flag to disable automatic elimination of structures containing rare earth elements.'
    )
    parser.add_argument(
        '--keep_equivalent',
        action='store_true',
        help='A flag to disable automatic structure matching and elimination of duplicates.'
    )

    args: Namespace = parser.parse_args()

    if args.output == '[filename]_out.cif':
        args.output = args.filename.replace('.cif', '_out.cif')

    assert_args(args)


    # MAIN BLOCK

    from screening_pipeline.utils import read_cif, write_cif
    from screening_pipeline.utils import batch_symmetrizer
    from screening_pipeline.utils import remove_equivalent

    # Extraction des données CIF et conversion en structures

    structures, nbr_rare_gas_structs, nbr_rare_earth_structs = read_cif(
        filename=args.filename,
        workers=args.workers,
        keep_rare_gases=args.keep_rare_gases, 
        keep_rare_earths=args.keep_rare_earths
    )

    nbr_loaded_structs = len(structures)
    nbr_total_structs  = nbr_loaded_structs + nbr_rare_gas_structs + nbr_rare_earth_structs
    assert nbr_loaded_structs > 0, 'No structure could be parsed from given data'

    print(f'{nbr_total_structs} structures detected in total')

    if not args.keep_rare_gases:
        print(f'{nbr_rare_gas_structs} structures containing rare gases were ignored')
    
    if not args.keep_rare_earths:
        print(f'{nbr_rare_earth_structs} structures containing rare earths were ignored')

    print(f'{nbr_loaded_structs} structures are kept for further processing')
    
    # Calcul de la symétrie d'espace des structures

    symmetrized_structs = list(filter(
        None, 
        batch_symmetrizer(
            structures=structures, 
            valid_tol=args.valid_tol, 
            symprec=args.precision, 
            angle_tolerance=args.angleprec, 
            workers=args.workers
        )
    ))

    print(f'{len(symmetrized_structs)} structures were symmetrized')

    if args.valid_tol > 0.0:
        nbr_not_valid = nbr_loaded_structs - len(symmetrized_structs)
        print(f'{nbr_not_valid} structures having too close atoms were discarded')

    # Comparaison des structures pour éliminer les doublons

    kept_structs, nbr_equivalent = remove_equivalent(
            structures=symmetrized_structs, 
            workers=args.workers, 
            keep_equivalent=args.keep_equivalent
    )

    nbr_unique_structs = len(kept_structs)

    if not args.keep_equivalent:
        print(f'{nbr_unique_structs} unique structures detected')
        print(f'{nbr_equivalent} duplicates were discarded')

    # Ecriture du fichier CIF symétrisé et épuré des structures indésirables
    write_cif(
        filename=args.output,
        structures=kept_structs,
        workers=args.workers,
    )

    # Calcul du temps total pris par la procédure

    stop = datetime.now()

    print(' ')
    print('------------------------------')
    print(' ')
    print('SUMMARY OF THE CALCULATION')
    print(' ')
    print(f'{nbr_total_structs} structures detected in total, including:')
    print(f'- {nbr_unique_structs} unique structures')

    if not args.keep_rare_gases:
        print(f'- {nbr_rare_gas_structs} structures containing rare gases')

    if not args.keep_rare_earths:
        print(f'- {nbr_rare_earth_structs} structures containing rare earths')
    
    if args.valid_tol > 0.0:
        print(f'{nbr_not_valid} sttructures with atoms that are too close (< {args.valid_tol} Angströms)')
    
    if not args.keep_equivalent:
        print(f'- {nbr_equivalent} structures that are duplicates')

    print(' ')
    print(f"Output results written in '{args.output}'")
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()
