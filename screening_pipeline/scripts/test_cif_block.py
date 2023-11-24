from argparse import ArgumentParser
from pymatgen.io.cif import CifWriter
from screening_pipeline.utils import read_cif


def main():

    # ARGUMENTS PARSING BLOCK

    prog_name = 'test_cif_block.py'
    prog_desc = '''
        testing bypassing techniques to avoid recalculation of symmetry on SymmetrizedStructures'''

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_desc
    )

    parser.add_argument(
        'filename', 
        type=str, 
        help='file to test.'
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

    args = parser.parse_args()


    # MAIN BLOCK

    structures = read_cif(
        filename=args.filename, 
        symprec=args.precision, 
        angle_tolerance=args.angleprec, 
        workers=args.workers, 
        keep_rare_gases=args.keep_rare_gases
    )
    structure  = structures[0]
    cif_writer = CifWriter(structure)
    cif_file   = cif_writer.cif_file
    cif_data   = cif_file.data
    cif_block  = cif_data['TiO2']
    block_data = cif_block.data
    loops_data = cif_block.loops
    print(f"block_data:\n{block_data}")
    print(f"loops_data:\n{loops_data}")
    print(type(block_data), type(loops_data))
    print(block_data['_symmetry_space_group_name_H-M'], type(block_data['_symmetry_space_group_name_H-M']))
    print(block_data['_symmetry_Int_Tables_number'], type(block_data['_symmetry_Int_Tables_number']))
    block_data['_symmetry_space_group_name_H-M'] = 'Pmmm'
    block_data['_symmetry_Int_Tables_number'] = 47
    print(block_data['_symmetry_space_group_name_H-M'], type(block_data['_symmetry_space_group_name_H-M']))
    print(block_data['_symmetry_Int_Tables_number'], type(block_data['_symmetry_Int_Tables_number']))

'''
def struct_to_cif_str(structure: Structure | SymmetrizedStructure):

    if not isinstance(structure, Structure):
        raise TypeError('Cannot write CIF data for a non-structure object.')

    cif_writer = CifWriter(structure)

    if not isinstance(structure, SymmetrizedStructure):
        cif_str = '# symmetrize.py: unable to find symmetry\n' + str(cif_writer)
        return cif_str
    
    cif_block_list = list(cif_writer.cif_file.data.values())
    cif_block      = cif_block_list[0]
    data_dict      = cif_block.data
    data_dict['_symmetry_space_group_name_H-M'] = structure.get_space_group_symbol()
    data_dict['_symmetry_Int_Tables_number'] = structure.get_space_group_number()
    symm_ops = structure.get_symmetry_operations()
    str_ops  = [op.as_xyz_string() for op in symm_ops]
    data_dict['_symmetry_equiv_pos_site_id'] = [f'{i}' for i in range(1, len(str_ops) + 1)]
    data_dict['_symmetry_equiv_pos_as_xyz'] = str_ops
    cif_str = str(cif_block)
    return cif_str
'''
if __name__=='__main__':
    main()