#!/usr/bin/python
'''
A script that calculates symmetry spacegroup for structures in a CIF file using pymatgen.
'''

##################################################
# SYSTEM I/O MODULES

from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, Namespace


def assert_args(args: Namespace):
    '''
    Input arguments verification.

    Args:
        args (NamedTuple): namespace of the parsed arguments.
    '''

    print(' - INPUT ARGUMENTS - ')
    print(f'filename: {args.filename}')
    print(f'output: {args.output}')
    print(f'equivalent: {args.equivalent}')

    assert (
        args.filename.endswith('.cif')
        and args.output.endswith('.cif')
        and (args.equivalent is None or args.equivalent.endswith('.cif'))
    ), 'some arguments formats are not supported, please only use CIF format'

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
    print('----------------------------------------')
    print(' ')


def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = 'symmetrize'
    prog_description = 'A script that calculates symmetry spacegroup for structures in a CIF file using pymatgen.'
    prog_missing_steps = ''''''
    helper_format = ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_description, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        'filename',
        type=str,
        help='name of the file containing the input structures'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default='[filename]_out.cif',
        help='name of output file containing the unique structures',
    )
    parser.add_argument(
        '-e',
        '--equivalent',
        type=str,
        default=None,
        help='name of the output file containing the duplicated structures',
    )
    parser.add_argument(
        '-p',
        '--precision',
        type=float,
        default=0.01,
        help='Fractional coordinates tolerance for symmetry finding',
        metavar='float',
    )
    parser.add_argument(
        '-a',
        '--angleprec',
        type=float,
        default=5.0,
        help='Angle tolerance for symmetry finding in degrees',
        metavar='float',
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create',
        metavar='int',
    )
    parser.add_argument(
        '--keep_rare_gases',
        action='store_true',
        help='Pass this flag to disable automatic elimination of structures containing rare gases'
    )
    parser.add_argument(
        '--keep_equivalent',
        action='store_true',
        help='Pass this flag to disable automatic structure matching and elimination of duplicates'
    )

    args: Namespace = parser.parse_args()

    if args.output == '[filename]_out.cif':
        args.output = args.filename.replace('.cif', '_out.cif')

    assert_args(args)


    # MAIN BLOCK

    from screening_pipeline.utils import read_cif, write_cif, remove_equivalent

    # Extraction des structures sous forme de strings

    symmetrized_structs = read_cif(
        args.filename,
        symprec=args.precision,
        angle_tolerance=args.angleprec,
        workers=args.workers,
        keep_rare_gases=args.keep_rare_gases,
    )

    assert len(symmetrized_structs) > 0, 'No structure could be parsed from given data'

    print(f'{len(symmetrized_structs)} structures loaded')

    # Comparaison des structures pour éliminer les doublons

    if args.keep_equivalent:
        kept_structs = symmetrized_structs
        duplicated_struct = []
    else:
        kept_structs, duplicated_struct = remove_equivalent(
            structures=symmetrized_structs, 
            workers=args.workers, 
            keep_equivalent=False
        )

        print(f'{len(kept_structs)} unique structures detected')

    #  Recalcul des symétries avec PyMatGen (pour prise en compte par CifWriter) et écriture du fichier de sortie

    write_cif(
        args.output,
        kept_structs,
        symprec=args.precision,
        angle_tolerance=args.angleprec,
        workers=args.workers,
    )

    if args.equivalent is not None:
        write_cif(
            args.equivalent,
            duplicated_struct,
            symprec=args.precision,
            angle_tolerance=args.angleprec,
            workers=args.workers,
        )

    # Calcul du temps total pris par la procédure

    stop = datetime.now()
    count_unique = len(kept_structs)
    count_duplicated = len(duplicated_struct)
    print(' ')
    print('------------------------------')
    print(' ')
    print('SUMMARY OF THE CALCULATION')
    print(' ')
    print(f'{count_duplicated+count_unique} structures detected in total, including:')
    print(f'- {count_duplicated} duplicated structures')
    print(f'- {count_unique} unique structures')
    print(' ')
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()
