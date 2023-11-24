"""
Functions to load and write cif structure with multiple workers.
"""


##################################################
# SYSTEM I/O MODULES

from typing import Tuple, List, Optional
from contextlib import redirect_stdout, redirect_stderr
import re

##################################################
# OPTIMIZATION MODULES

from itertools import filterfalse, count
from tqdm.contrib.concurrent import process_map

##################################################
# PYTHON MATERIALS GENOMICS MODULE

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from screening_pipeline.utils import has_rare_gas, discard_rare_gas_structures
from screening_pipeline.utils.redirect import redirect_c_stdout, redirect_c_stderr

##################################################


def retry_get_symmetrized_structure(
        structure: Structure,
        symprec: Optional[float] = None,
        angle_tolerance: float = 5.0,
    ) -> Structure:
    
    for precision_factor in [2, 3, 5, 10]:

        symmetrizer = SpacegroupAnalyzer(
            structure=structure, 
            symprec=precision_factor*symprec, 
            angle_tolerance=precision_factor*angle_tolerance
        )

        try:
            sym_struct = symmetrizer.get_symmetrized_structure()

        except TypeError:

            if symmetrizer.get_symmetry_dataset() is not None:
                return structure
            
        else:
            return sym_struct
    
    return structure

def extract_cif_from_file(filename: str) -> List[str]:
    """
    Separates concatenated structures from a cif file.

    Args:
        filename (str): Path to a cif file.
        keep_rare_gases (bool): Whether the structures containing rare gases should be kept or not. Default is false.
    Returns
        List[str]: List of CIF strings, each one representing a single structure.
    """

    # load file
    with open(filename, "r") as input_file:
        data = input_file.read()

    # match structures data with regular expression
    matcher = re.compile(r"^data.*?$(?=\ndata|\Z)", re.MULTILINE | re.DOTALL)
    full_structs_list = matcher.findall(data)

    return full_structs_list


def cif_str_to_struct(
    cif_str: str,
    symprec: Optional[float] = None,
    angle_tolerance: float = 5.0,
) -> Structure:
    """
    Parses data from a cif formatted string and converts it to a structure object from pymatgen. This function uses a spacegroup analyser from pymatgen and spglib in the backend.

    Args:
        cif_str (str): The string of a structure encoded in the cif format.
        symprec (float): Distance tolerance for symmetry search.
        angle_tolerance (float): Angle tolerance for symmetry search.
    Returns
        A Pymatgen Structure object if symprec is None, a Pymatgen SymmetrizedStructure object if not.
    """

    with redirect_c_stdout(None), redirect_c_stderr(None):
        parsed_str = CifParser.from_str(cif_string=cif_str)
        struct_list = parsed_str.get_structures()
        struct = struct_list[0]

        if symprec is None:
            return struct

        sym_struct = SpacegroupAnalyzer(
            struct, symprec, angle_tolerance
            )

        try:
            sym_struct = sym_struct.get_symmetrized_structure()

        except TypeError:

            if sym_struct.get_symmetry_dataset() is None:
                
                return retry_get_symmetrized_structure(
                    structure=struct, 
                    symprec=symprec, 
                    angle_tolerance=angle_tolerance
                )

            return struct

        return sym_struct

def _cif_str_to_struct_fn(args):
    return cif_str_to_struct(*args)


def struct_to_cif_str(
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

    return cif_str


def _struct_to_cif_str_fn(args):
    return struct_to_cif_str(*args)


def read_cif(
    filename: str,
    symprec: Optional[float] = None,
    angle_tolerance: float = 5.0,
    workers: int = 1,
    keep_rare_gases: bool = False,
) -> List[Structure]:
    """
    Read multiple structures from a cif file and decode them using multiprocess.

    Args:
        filename (str): Name of the input file.
        symprec (float): Distance tolerance for symmetry search.
        angle_tolerance (float): Angle tolerance for symmetry search.
        workers (int): Number of workers used.
        keep_rare_gases (bool): Whether the structures containing rare gases should be kept or not. Default is false.
    Returns
        List[Structure]: Returns the structures in a list.
    """

    struct_strings = extract_cif_from_file(filename)

    if not keep_rare_gases:
        struct_strings = discard_rare_gas_structures(struct_strings)

    nbr_struct = len(struct_strings)
    assert nbr_struct > 0, "No structure data found in provided file"

    # Obtention des structures PyMatGen à partir des données et calcul de la symétrie

    chunksize = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    def feed_args(structures, symprec, angle_tolerance) -> List[Tuple]:
        return [(struct, symprec, angle_tolerance) for struct in structures]

    return process_map(
        _cif_str_to_struct_fn,
        feed_args(struct_strings, symprec, angle_tolerance),
        max_workers=workers,
        chunksize=chunksize,
        desc="load and search symmetries",
    )


def write_cif(
    filename: str,
    structures: List[Structure],
    symprec: Optional[float] = None,
    angle_tolerance: float = 5.0,
    workers: int = 1,
):
    """
    Write multiple structures to a cif file and encode them using multiprocess.

    Args:
        filename (str): Name of the input file.
        structures (List[Structure]): The structures to encode.
        symprec (float): Distance tolerance for symmetry search.
        angle_tolerance (float): Angle tolerance for symmetry search.
        workers (int): Number of workers used.
    """
    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    def feed_args(structures, symprec, angle_tolerance) -> List[Tuple]:
        return [(struct, symprec, angle_tolerance) for struct in structures]

    encoded_cif = process_map(
        _struct_to_cif_str_fn,
        feed_args(structures, symprec, angle_tolerance),
        max_workers=workers,
        chunksize=chunksize,
        desc="convert to cif format",
    )

    with open(filename, "wt") as out_file:
        out_file.write("\n".join(encoded_cif))
