#!/usr/bin/python

'''
A script to determine material fundamental band gap from VASP relaxations and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (reference 32 in screening_pipeline/Bibliography)
'''


from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

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
        help='The file conaining structure data to calculate band gap from. CIF format only.'
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
    print(f"elapsed time: {stop-start}")


if __name__ == '__main__':
    main()