#!/usr/bin/python

'''
A script to determine material valence and conduction band edge positions.

Reference of the method:
    Y. Wu, M.K.Y. Chan, and G. Ceder, Phys. Rev. B, 83, 235301 (2011)
    (reference 27 in screening_pipeline/Bibliography)
'''


from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, Namespace


def assert_args(args: Namespace):
    '''
    Input arguments verification.

    Args:
        args (NamedTuple): namespace of the parsed arguments.
    '''

    assert args.filename.endswith('.cif'), 'Input file must be in CIF format'
    assert args.workers >= 1, 'the number of workers cannot be negative or zero'

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK
    
    prog_name = 'band_edge_pos_screening'
    prog_desc = '''
        A script to determine material valence and conduction band edge positions.

        Reference of the method:
            Y. Wu, M.K.Y. Chan, and G. Ceder, Phys. Rev. B, 83, 235301 (2011)
            (reference 27 in screening_pipeline/Bibliography)
        '''
    prog_missing_steps = '''
        Missing steps to complete this script:
            - Aknowledge the steps of the method
            - Add in-line command arguments
        '''
    helper_format = ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_desc, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        'filename', 
        type=str, 
        help='The file conaining structure data to calculate band edge positions from. CIF format only.'
    )
    parser.add_argument(
        '-p', '--path',
        type=str,
        default='./Vasp_input_sets/',
        help='Directory to write VASP input files in (created if it does not exist).',
        metavar='str'
    )
    parser.add_argument(
        '-w', '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create',
        metavar='int'
    )

    args: Namespace = parser.parse_args()

    assert_args(args)


    # MAIN BLOCK

    from screening_pipeline.utils import read_cif

    structures = read_cif(
        filename=args.filename, 
        workers=args.workers, 
        keep_rare_gases=args.keep_rare_gases
    )
    
    nbr_struct: int   = len(structures)
    chunksize: int    = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
    generic_path: str = args.path

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()