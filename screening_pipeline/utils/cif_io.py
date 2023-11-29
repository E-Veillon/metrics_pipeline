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
from functools import partial

########################################
# PYTHON MATERIALS GENOMICS PACKAGE

from pymatgen.core.structure import Structure, PeriodicSite
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
        cif_str = '# symmetrize: unable to find symmetry\n' + str(cif_writer)
        return cif_str
    
    # Extract the data dict from CifWriter
    cif_block_list      = list(cif_writer.cif_file.data.values())
    cif_block           = cif_block_list[0]
    data_dict           = cif_block.data
    # Extract symmetry infos from the SymmetrizedStructure object
    struct_spg          = structure.spacegroup
    xyz_ops             = [op.as_xyz_string() for op in struct_spg]
    equiv_sites         = structure.equivalent_sites
    # Get the ordering for positions matrix
    partial_sort        = partial(sorted, key=lambda s: tuple(abs(x) for x in s.frac_coords))
    unique_sites        = [(partial_sort(sites)[0],len(sites)) for sites in equiv_sites]
    partial_sort_2      = partial(sorted, key=lambda t: (t[0].species.average_electroneg,-t[1],t[0].a,t[0].b,t[0].c))
    sorted_unique_sites = partial_sort_2(unique_sites)


    data_dict['_symmetry_space_group_name_H-M'] = struct_spg.int_symbol
    data_dict['_symmetry_Int_Tables_number']    = struct_spg.int_number
    data_dict['_symmetry_equiv_pos_site_id']    = [f'{idx}' for idx in range(1, len(xyz_ops) + 1)]
    data_dict['_symmetry_equiv_pos_as_xyz']     = xyz_ops
    #TODO: prendre en compte les équivalences entre sites pour réduire la matrice des coordonnées
    for site, mult in sorted_unique_sites:
        for specie, occupancy in site.species.items():
            atom_site_type_symbol.append(str(specie))
            atom_site_symmetry_multiplicity.append(f"{mult}")
            atom_site_fract_x.append(format_str.format(site.a))
            atom_site_fract_y.append(format_str.format(site.b))
            atom_site_fract_z.append(format_str.format(site.c))
            site_label = site.label if site.label != site.species_string else f"{specie.symbol}{count}"
            atom_site_label.append(site_label)
            atom_site_occupancy.append(str(occupancy))
            count += 1
    # Ajout des variables au dictionnaire de données
    block["_atom_site_type_symbol"] = atom_site_type_symbol
    block["_atom_site_label"] = atom_site_label
    block["_atom_site_symmetry_multiplicity"] = atom_site_symmetry_multiplicity
    block["_atom_site_fract_x"] = atom_site_fract_x
    block["_atom_site_fract_y"] = atom_site_fract_y
    block["_atom_site_fract_z"] = atom_site_fract_z
    block["_atom_site_occupancy"] = atom_site_occupancy

    cif_str = str(cif_block)

    return cif_str

def read_cif(
    filename: str,
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
    workers: int = 1,
) -> None:
    """
    Encode multiple structures in CIF format and write them in a file using multiprocess.

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
            desc="Writing data in cif format"
        ))

    with open(filename, "wt") as out_file:
        out_file.write("\n".join(encoded_cif))
