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

from itertools import filterfalse
from tqdm.contrib.concurrent import process_map

##################################################
# PYTHON MATERIALS GENOMICS MODULE

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from screening_pipeline.utils import has_rare_gas
from screening_pipeline.utils.redirect import redirect_c_stdout, redirect_c_stderr

##################################################


def extract_cif_from_file(
    filename: str, keep_rare_gases: bool = False
) -> Tuple[List[str], List[str]]:
    """
    Separates concatenated structures from a cif file.

    Args:
        filename (str): Path to a cif file.
        keep_rare_gases (bool): Whether the structures containing rare gases should be kept or not. Default is false.
    Returns
        Tuple[List[str], List[str]]: Returns the list of filtered structures and the list of removed structures.
    """

    # load file
    with open(filename, "r") as input_file:
        data = input_file.read()

    # match structures data with regular expression
    matcher = re.compile(r"^data.*?$(?=\ndata|\Z)", re.MULTILINE | re.DOTALL)
    full_structs_list = matcher.findall(data)

    if keep_rare_gases:
        return full_structs_list, []

    # filter structures if they contain rare gases
    rare_gas_structures = list(filter(has_rare_gas, full_structs_list))
    full_structs_list = list(filterfalse(has_rare_gas, full_structs_list))

    return full_structs_list, rare_gas_structures


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
            ).get_symmetrized_structure()

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
        cif_str = str(
            CifWriter(struct=struct, symprec=symprec, angle_tolerance=angle_tolerance)
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

    structures, _ = extract_cif_from_file(
        filename, keep_rare_gases
    )
    nbr_struct = len(structures)
    assert nbr_struct > 0, "No structure data found in provided file"

    # Obtention des structures PyMatGen à partir des données et calcul de la symétrie

    chunksize = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    def feed_args(structures, symprec, angle_tolerance):
        return [(struct, symprec, angle_tolerance) for struct in structures]

    return process_map(
        _cif_str_to_struct_fn,
        feed_args(structures, symprec, angle_tolerance),
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

    def feed_args(structures, symprec, angle_tolerance):
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
