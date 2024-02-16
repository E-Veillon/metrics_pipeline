#!/usr/bin/python
'''
A script that calculates symmetry spacegroup for structures in a CIF file using pymatgen.
'''

########################################
# SYSTEM I/O MODULES

import os
from datetime import datetime
from argparse import ArgumentParser, RawTextHelpFormatter, Namespace

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace):
    '''
    Apply viability assertions on input arguments and prints a summary of them.

    Parameters:
        args (Namespace): namespace of the parsed arguments.
    '''

    print(' - I/O ARGUMENTS - ')
    print(f'INPUT FILE: {args.input_file}')
    print(f'OUTPUT FILE: {args.output}')

    assert os.path.isfile(args.input_file) and args.input_file.endswith('.cif'), \
    f'{args.input_file} is not an existing CIF file.'

    assert args.output.endswith('.cif'), \
    f'Output file must be of CIF format.'

    print(' ')
    print('------------------------------')
    print(' ')
    print(' - ACTIVATED FEATURES - ')
    print(f'CHECK RARE GASES: {not args.no_rare_gas_check}')
    print(f'CHECK RARE EARTHS: {not args.no_rare_earth_check}')
    print(f'CHECK INTERATOMIC DISTANCES: {not args.no_dist_check}')
    print(f'* Distance tolerance: {args.dist_tolerance} Angstroms {"(ignored)" if args.no_dist_check else ""}')

    if not args.no_dist_check:
        assert args.dist_tolerance > 0.0, \
        'Interatomic distance tolerance must be positive.'

    print(f'SYMMETRIZATION: {not args.no_symmetrization}')
    print(f'* Fractional coordinates tolerance: {args.symprec} {"(ignored)" if args.no_symmetrization else ""}')

    if not args.no_symmetrization:
        assert (0.0 <= args.symprec <= 0.5), \
        'Fractional coordinates tolerance must be between 0.0 and 0.5 to retain some reliability.'

    print(f'* Angles tolerance: {args.angleprec} degrees {"(ignored)" if args.no_symmetrization else ""}')

    if not args.no_symmetrization:
        assert (0.0 <= args.angleprec <= 30.0), \
        'Angles tolerance must be between 0.0 and 20.0 degrees to retain some reliability.'

    print(f'STRUCTURE MATCHING: {not args.no_equiv_match}')
    print(f'NUMBER OF WORKERS: {args.workers}')

    assert args.workers >= 1, 'The number of workers must be positive.'

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
    #prog_missing_steps = '''
    #    Missing steps to complete this script: None
    #    '''
    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_description, 
        #epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        'input-file',
        type=str,
        help='Name or path to the CIF file containing structure data to process.'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default='[input_file]_out.cif',
        help='''Name or path to wanted CIF output file where processed structures will be stored 
                (default: %(default)s).''',
    )
    parser.add_argument(
        '--no-rare-gas-check',
        action='store_true',
        help='A flag to disable elimination of structures containing rare gas elements.'
    )
    parser.add_argument(
        '--no-rare-earth-check',
        action='store_true',
        help='A flag to disable elimination of structures containing f-block elements.'
    )
    parser.add_argument(
        '--no-dist-check', 
        action='store_true', 
        help='A flag to disable structures interatomic distances checking.'
    )
    parser.add_argument(
        '-d', 
        '--dist-tolerance', 
        type=float, 
        default=0.5, 
        help='''Tolerance for checking interatomic distances in Angstroms.
                Structures containing atoms that are closer than this value will be discarded.
                (Default: %(default)s Angstroms).''', 
        metavar='float', 
        dest='valid_tol'
    )
    parser.add_argument(
        '--no-symmetrization', 
        action='store_true', 
        help='A flag to disable search of structures symmetry space groups.'
    )
    parser.add_argument(
        '-s',
        '--symprec',
        type=float,
        default=0.01,
        help='Fractional coordinates tolerance for symmetry finding (Default: %(default)s).',
        metavar='float',
    )
    parser.add_argument(
        '-a',
        '--angleprec',
        type=float,
        default=5.0,
        help='Angle tolerance for symmetry finding in degrees (Default: %(default)s degrees).',
        metavar='float',
    )
    parser.add_argument(
        '--no-equiv-match',
        action='store_true',
        help='A flag to disable structure matching and elimination of duplicates.'
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to spawn to speed up the process (Default: %(default)s).',
        metavar='int',
    )

    args: Namespace = parser.parse_args()

    if args.output == '[filename]_out.cif':
        args.output = args.input_file.replace('.cif', '_out.cif')

    assert_args(args)

    input_file      = args.input_file
    output_file     = args.output
    keep_rare_gas   = args.no_rare_gas_check
    keep_rare_earth = args.no_rare_earth_check
    no_dist_check   = args.no_dist_check
    dist_tol        = args.dist_tolerance
    no_symmetry     = args.no_symmetrization
    symprec         = args.symprec
    angleprec       = args.angleprec
    no_equiv_match  = args.no_equiv_match
    workers         = args.workers


    # MAIN BLOCK

    from screening_pipeline.utils.cif_io import read_cif, write_cif
    from screening_pipeline.utils.data_process import check_interatomic_distances
    from screening_pipeline.utils.spacegroup import batch_symmetrizer
    from screening_pipeline.utils.matcher import remove_equivalent

    # Extraction des données CIF et conversion en structures

    structures, nbr_rare_gas_structs, nbr_rare_earth_structs = read_cif(
        filename=input_file,
        workers=workers,
        keep_rare_gases=keep_rare_gas, 
        keep_rare_earths=keep_rare_earth
    )

    nbr_loaded_structs = len(structures)
    nbr_total_structs  = nbr_loaded_structs + nbr_rare_gas_structs + nbr_rare_earth_structs
    assert nbr_loaded_structs > 0, 'No structure could be parsed from given data'

    print(f'{nbr_total_structs} structures detected in total')

    if not keep_rare_gas:
        print(f'{nbr_rare_gas_structs} structures containing rare gases were ignored')
    
    if not keep_rare_earth:
        print(f'{nbr_rare_earth_structs} structures containing rare earths were ignored')

    print(f'{nbr_loaded_structs} structures are kept for further processing')
    
    # Vérification des distances interatomiques

    if not no_dist_check:
        structures, nbr_not_valid = check_interatomic_distances(structures, valid_tol=dist_tol)
        print(f'{nbr_not_valid} structures having too close atoms were discarded')

    # Calcul de la symétrie d'espace des structures

    if no_symmetry: symmetrized_structs = structures
    else:
        symmetrized_structs = list(filter(
            None, 
            batch_symmetrizer(
                structures=structures,  
                symprec=args.precision, 
                angle_tolerance=args.angleprec, 
                workers=workers
            )
        ))

        print(f'{len(symmetrized_structs)} structures were symmetrized')

    # Comparaison des structures pour éliminer les doublons

    kept_structs, nbr_equivalent = remove_equivalent(
            structures=symmetrized_structs, 
            workers=workers, 
            keep_equivalent=no_equiv_match
    )

    nbr_unique_structs = len(kept_structs)

    if not no_equiv_match:
        print(f'{nbr_unique_structs} unique structures detected')
        print(f'{nbr_equivalent} duplicates were discarded')

    # Ecriture du fichier CIF symétrisé et épuré des structures indésirables
    write_cif(
        filename=output_file,
        structures=kept_structs,
        workers=workers,
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

    if not keep_rare_gas:
        print(f'- {nbr_rare_gas_structs} structures containing rare gases')

    if not keep_rare_earth:
        print(f'- {nbr_rare_earth_structs} structures containing rare earths')
    
    if not no_dist_check:
        print(f'- {nbr_not_valid} structures with too small interatomic distances')
    
    if not no_equiv_match:
        print(f'- {nbr_equivalent} structures that are duplicates')

    print(' ')
    print(f"Output results written in '{output_file}'")
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()
