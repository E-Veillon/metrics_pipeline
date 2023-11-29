'''
Functions to load and write CIF formatted data with multiple processes.
'''


########################################
# SYSTEM I/O MODULES

from typing import Tuple, List, Union
#from contextlib import redirect_stdout, redirect_stderr
import re

########################################
# OPTIMIZATION MODULES

from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIALS GENOMICS PACKAGE

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SymmetrizedStructure

########################################
# LOCAL MODULES

from screening_pipeline.utils import discard_rare_gas_structures
from screening_pipeline.utils.redirect import redirect_c_stdout, redirect_c_stderr

##################################################


def extract_cif_from_file(filename: str) -> List[str]:
    """
    Separates concatenated structures from a cif file.

    Parameters:
        filename (str):         Path to a cif file.

        keep_rare_gases (bool): Whether the structures containing rare gases should be kept or not. 
                                Defaults to false.
    
    Returns:
        List[str]: List of CIF strings, each one representing a single structure.
    """

    # load file
    with open(filename, "r") as input_file:
        data = input_file.read()

    # match structures data with regular expression
    matcher = re.compile(r"^data.*?$(?=\ndata|\Z)", re.MULTILINE | re.DOTALL)
    cif_str_structs = matcher.findall(data)

    return cif_str_structs


def cif_str_to_struct(cif_str: str) -> Structure:
    """
    Parses data from a cif formatted string and converts it to a pymatgen Structure object.

    Parameters:
        cif_str (str): The string of a structure encoded in the cif format.

    Returns:
        A Pymatgen Structure object.
    """

    with redirect_c_stdout(None), redirect_c_stderr(None):
        parsed_str  = CifParser.from_str(cif_string=cif_str)
        struct_list = parsed_str.parse_structures()
        structure   = struct_list[0]
        return structure

'''def struct_to_cif_str(
    struct: Structure,
    symprec: Optional[float] = None,
    angle_tolerance: float = 5.0,
) -> str:
    """
    Calculates Hermann-Mauguin's spacegroup in a pymatgen's Structure object and converts it into cif formatted data string.

    Args:
        struct (Structure): The structure to encode.
        symprec (float): Distance tolerance for symmetry search.
        angle_tolerance (float): Angle tolerance for symmetry search.
    Returns
        str: Returns the encoded structure.
    """

    with redirect_c_stdout(None), redirect_c_stderr(None):

        try:
            cif_str = str(
                CifWriter(struct=struct, symprec=symprec, angle_tolerance=angle_tolerance)
            )
        except TypeError:
            cif_str = str(
                CifWriter(struct=struct)
            )

    return cif_str'''

def struct_to_cif_str(
        structure: Union[Structure, SymmetrizedStructure]
    ) -> str:
    '''
    Converts a pymatgen Structure object into a CIF formatted string.

    Parameters:
        structure (Structure): Structure object to convert.
    
    Returns:
        str: CIF formatted string.
    '''

    if not isinstance(structure, Structure):
        raise TypeError('Cannot write CIF data for a non-structure object.')

    cif_writer = CifWriter(structure)

    if not isinstance(structure, SymmetrizedStructure):
        cif_str = '# symmetrize.py: unable to find symmetry\n' + str(cif_writer)
        return cif_str
    
    cif_block_list = list(cif_writer.cif_file.data.values())
    cif_block      = cif_block_list[0]
    data_dict      = cif_block.data
    symm_ops       = structure.get_symmetry_operations()
    str_ops        = [op.as_xyz_string() for op in symm_ops]

    data_dict['_symmetry_space_group_name_H-M'] = structure.get_space_group_symbol()
    data_dict['_symmetry_Int_Tables_number']    = structure.get_space_group_number()
    data_dict['_symmetry_equiv_pos_site_id']    = [f'{idx}' for idx in range(1, len(str_ops) + 1)]
    data_dict['_symmetry_equiv_pos_as_xyz']     = str_ops
    #TODO: prendre en compte les équivalences entre sites pour réduire la matrice des coordonnées
    cif_str = str(cif_block)

    return cif_str

def read_cif(
    filename: str,
    #symprec: Optional[float] = None,
    #angle_tolerance: float = 5.0,
    workers: int = 1,
    keep_rare_gases: bool = False,
) -> Tuple[List[Structure], int]:
    """
    Reads a cif file containing concatenated structures data and decode them using multiprocess.

    Parameters:
        filename (str):         Name of the input CIF file.
        
        workers (int):          Number of processes to use in parallel.
        
        keep_rare_gases (bool): Whether structures containing rare gases should be kept. 
                                Defaults to false.
    
    Returns:
        List[Structure]: Decoded Structure objects in a list.
        Int: Number of structures containing rare gases discarded. 
    """

    struct_strings = extract_cif_from_file(filename)

    if not keep_rare_gases:
        struct_strings, nbr_rare_gas_structs = discard_rare_gas_structures(struct_strings)

    nbr_struct = len(struct_strings)
    assert nbr_struct > 0, "No structure data found in provided file"

    chunksize = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    #def feed_args(structures, symprec, angle_tolerance) -> List[Tuple]:
    #    return [(struct, symprec, angle_tolerance) for struct in structures]

    structs_list = list(process_map(
        cif_str_to_struct,
        struct_strings,
        max_workers=workers,
        chunksize=chunksize,
        desc="load and read data",
    ))

    return structs_list, nbr_rare_gas_structs

def write_cif(
    filename: str,
    structures: List[Structure],
    #symprec: Optional[float] = None,
    #angle_tolerance: float = 5.0,
    workers: int = 1,
) -> None:
    """
    Encode multiple structures in CIF formatand write them in a file using multiprocess.

    Parameters:
        filename (str):               Name of the input file.

        structures (List[Structure]): The structures to encode.

        workers (int):                Number of parallel processes to use.
    """

    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    #def feed_args(structures, symprec, angle_tolerance) -> List[Tuple]:
    #    return [(struct, symprec, angle_tolerance) for struct in structures]

    encoded_cif = list(
        process_map(
            struct_to_cif_str, 
            structures, 
            #feed_args(structures, symprec, angle_tolerance),
            max_workers=workers, 
            chunksize=chunksize, 
            desc="converting to cif format"
        ))

    with open(filename, "wt") as out_file:
        out_file.write("\n".join(encoded_cif))
