#!/usr/bin/python

'''
A script to determine material valence and conduction band edge positions.

Reference of the method:
    Y. Wu, M.K.Y. Chan, and G. Ceder, Phys. Rev. B, 83, 235301 (2011)
    (reference 27 in screening_pipeline/Bibliography)
'''


from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

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

    # MAIN BLOCK

    

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()