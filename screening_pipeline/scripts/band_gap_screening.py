#!/usr/bin/python

'''
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (reference 32 in screening_pipeline/Bibliography)
'''


from datetime import datetime
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter


def assert_args(args: Namespace):
    '''
    Input arguments verification.

    Parameters:
        args (Namespace): namespace of the parsed arguments.
    '''

    assert args.filename.endswith('.cif'), 'Input file must be in CIF format'
    assert args.workers >= 1, 'the number of workers cannot be negative or zero'

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
    helper_format = RawDescriptionHelpFormatter

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
        help='Directory to write VASP input files in (created if it does not exist) (default = ./Vasp_input_sets/).',
        metavar='str'
    )
    parser.add_argument(
        '-w', '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create (default = 1).',
        metavar='int'
    )

    args: Namespace = parser.parse_args()

    assert_args(args)


    # MAIN BLOCK

    from screening_pipeline.utils.cif_io import read_cif
    from screening_pipeline.utils.periodic_table import get_all_valence_electrons, get_delta_sol_el_ratio
    from screening_pipeline.utils.vasp_io import vasp_input_files_settings, vasp_launcher

    structures = read_cif(
        filename=args.filename, 
        workers=args.workers, 
        keep_rare_gases=args.keep_rare_gases
    )

    nbr_struct: int   = len(structures)
    chunksize: int    = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
    generic_path: str = args.path

    for struct in structures:
        N_0 = get_all_valence_electrons(struct)
        n   = get_delta_sol_el_ratio(struct, 'PBE', 'BEST')
        vasp_input = vasp_input_files_settings(struct)
        vasp_launcher(vasp_input, generic_path) # Relaxation -> E(N0)
        # Changer la valeur de la densité de charge dans CHGCAR :
        # Densité de charge ponctuelle n(r) = densité électronique ponctuelle rhô(r) x charge élémentaire e
        # Densité de charge CHGCAR nc(r) = densité de charge ponctuelle n(r) x Vol. maille V(maille)
        # Ccl: nc(r) représente la densité de charge par maille, on peut donc enlever/ajouter
        # une densité de charge n x e directement à nc(r) pour faire les calculs statiques E(N0 - n) et E(N0 + n)

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()