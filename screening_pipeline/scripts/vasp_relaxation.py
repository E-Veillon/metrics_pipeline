#!/usr/bin/python
'''
A script that writes and performs VASP relaxation on structure data read from CIF files,
using pymatgen as a setting interface between raw data and VASP.
The relaxation results then may be used in other scripts for material properties analysis.
'''

########################################
# SYSTEM I/O MODULES

from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter

########################################
# OPTIMIZATION MODULES

from functools import partial
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIALS GENOMICS PACKAGE

from pymatgen.core.structure import Structure, SiteCollection

########################################
# LOCAL MODULES

from screening_pipeline.utils.cif_io import read_cif
from screening_pipeline.utils.vasp_io import vasp_relaxation_settings, vasp_launcher

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:
    assert args.filename.endswith('.cif'), \
    'Input structure data must be in CIF format.'

    allowed_presets = {(
        'MITRelaxSet', 
        'MPRelaxSet', 
        'MPScanRelaxSet', 
        'MPHSERelaxSet', 
        'MPMetalRelaxSet', 
        'MVLRelax52Set', 
        'MVLScanRelaxSet'
    )}

    assert args.method in allowed_presets, '''
    The relaxation method must be one of the following pymatgen relaxation presets:
    MITRelaxSet, MPRelaxSet, MPScanRelaxSet, MPHSERelaxSet.'''

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = 'vasp_relaxation.py'
    prog_description = '''
    A script that writes and performs VASP relaxation on structure data read from CIF files,
    using pymatgen as a setting interface between raw data and VASP.
    The relaxation results then may be used in other scripts for material properties analysis.
    '''
    prog_missing_steps = '''
        Missing steps to complete this script: 
            - Read cif, transform to structures - OK
            - Setup VASP calculation
            - Write VASP input files
            - Run VASP on written directories
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
        help='The CIF file containing structure data to read.', 
        metavar='file.cif'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default='./',
        help='''Output directory where VASP files will be written.
                A subdirectory will be created in output directory
                for each structure found in file.cif.''', 
        metavar='outdir'
    )
    parser.add_argument(
        '-m', 
        '--method', 
        type=str, 
        default='MITRelaxSet', 
        help='''The pymatgen preset to use for VASP relaxation.
                More info on possible presets in pymatgen documentation:
                https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.''', 
        metavar='RelaxSet'
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create',
        metavar='int',
    )
    args: Namespace = parser.parse_args()

    assert_args(args)
    
    filename = args.filename

    if not args.output.endswith('/'):
        args.output = args.output + '/'

    outdir  = args.output
    method  = args.method
    workers = args.workers


    # MAIN BLOCK

    structures, _, _ = read_cif(filename=filename)
    
    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
    vasp_setup = partial(
        vasp_relaxation_settings, 
        method=method, 
    )
    
    vasp_inputs = list(process_map(
        vasp_setup, 
        structures, 
        workers=workers, 
        chunksize=chunksize
    ))

    stop = datetime.now()
    print(f'Elapsed time: {stop-start}')


if __name__=='__main__':
    main()