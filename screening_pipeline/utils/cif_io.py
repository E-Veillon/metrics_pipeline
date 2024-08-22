#!/usr/bin/python
"""
Functions to load and write CIF formatted data with multiple processes.
"""


import re
from typing import Tuple, List
from tqdm.contrib.concurrent import process_map

# PYTHON MATERIALS GENOMICS
from pymatgen.core import Structure, PeriodicSite
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SymmetrizedStructure

# LOCAL IMPORTS
from .common_asserts import check_type, check_num_value
from .file_io import check_file_format, check_file_or_dir
from .periodic_table import (
    discard_rare_gas_structures, discard_rare_earth_structures
)
from .redirect import redirect_c_stdout, redirect_c_stderr


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
    check_file_or_dir(filename, "file", allowed_formats="cif")

    # load file
    with open(filename, "rt", encoding="utf-8") as input_file:
        data = input_file.read()

    # match structures data with regular expression
    matcher = re.compile(r"^data.*?$(?=\ndata|\Z)", re.MULTILINE | re.DOTALL)
    cif_str_structs = matcher.findall(data)

    return cif_str_structs


########################################


def cif_str_to_struct(cif_str: str) -> Structure:
    """
    Parses data from a cif formatted string and converts it to a pymatgen Structure object.

    Parameters:
        cif_str (str): The string of a structure encoded in the cif format.

    Returns:
        A Pymatgen Structure object.
    """
    check_type(cif_str, "cif_str", (str,))

    with redirect_c_stdout(None), redirect_c_stderr(None):
        parsed_str  = CifParser.from_str(cif_string=cif_str)
        struct_list = parsed_str.parse_structures(primitive=False)
        structure   = struct_list[0]
        return structure


########################################


def struct_to_cif_str(
        structure: Structure,
        significant_figures: int = 8
    ) -> str:
    """
    Converts a pymatgen Structure object into a CIF formatted string.

    Parameters:
        structure (Structure): Structure object to convert.

        significant_figures (int): Number of decimal places to keep.
    
    Returns:
        str: CIF formatted string.
    """
    check_type(structure, "structure", (Structure,))
    check_type(significant_figures, "significant_figures", (int,))
    check_num_value(significant_figures, "significant_figures", ">=", 0)

    cif_writer = CifWriter(
        struct=structure,
        significant_figures=significant_figures
    )

    if not isinstance(structure, SymmetrizedStructure):
        cif_str = str(cif_writer)
        return cif_str

    # Extract the data dict from CifWriter
    cif_block_list = list(cif_writer.cif_file.data.values())
    cif_block = cif_block_list[0]
    data_dict = cif_block.data

    # Extract symmetry infos from the SymmetrizedStructure object
    xyz_ops = [op.as_xyz_str() for op in structure.spacegroup]
    equiv_sites = structure.equivalent_sites

    # Get the ordering for species listing

    # Inside a list of equivalent sites, we only keep the one with minimal
    # fractional coordinates, and take the number of equivalent sites as the multiplicity
    def _frac_coord_sort_key(s: SymmetrizedStructure):
        return tuple(abs(x) for x in s.frac_coords)

    unique_sites: List[Tuple[PeriodicSite, int]] = [
        (
            sorted(sites, key=_frac_coord_sort_key)[0],
            len(sites)
        ) for sites in equiv_sites
    ]

    # Between non-equivalent sites, we sort them firstly by ascending electronegativity,
    # then by descending multiplicity, then by ascending frac coordinates
    def _electroneg_sort_key(t: Tuple):
        return (t[0].species.average_electroneg, -t[1], t[0].a, t[0].b, t[0].c)

    sorted_unique_sites: List[Tuple[PeriodicSite, int]] = sorted(
        unique_sites, key=_electroneg_sort_key
    )

    atom_site_type_symbol   = []
    atom_site_symmetry_mult = []
    atom_site_fract_x       = []
    atom_site_fract_y       = []
    atom_site_fract_z       = []
    atom_site_label         = []
    atom_site_occupancy     = []

    format_str = f"{{:.{significant_figures}f}}"
    count = 0

    for site, mult in sorted_unique_sites:
        for specie, occupancy in site.species.items():
            atom_site_type_symbol.append(str(specie))
            atom_site_symmetry_mult.append(f"{mult}")
            atom_site_fract_x.append(format_str.format(site.a))
            atom_site_fract_y.append(format_str.format(site.b))
            atom_site_fract_z.append(format_str.format(site.c))
            site_label = (
                site.label if site.label != site.species_string
                else f"{specie.symbol}{count}"
            )
            atom_site_label.append(site_label)
            atom_site_occupancy.append(str(occupancy))
            count += 1

    data_dict["_symmetry_space_group_name_H-M"] = structure.spacegroup.int_symbol
    data_dict["_symmetry_Int_Tables_number"] = structure.spacegroup.int_number
    data_dict["_symmetry_equiv_pos_site_id"] = [
        f"{idx}" for idx in range(1, len(xyz_ops) + 1)
    ]
    data_dict["_symmetry_equiv_pos_as_xyz"] = xyz_ops
    data_dict["_atom_site_type_symbol"] = atom_site_type_symbol
    data_dict["_atom_site_label"] = atom_site_label
    data_dict["_atom_site_symmetry_multiplicity"] = atom_site_symmetry_mult
    data_dict["_atom_site_fract_x"] = atom_site_fract_x
    data_dict["_atom_site_fract_y"] = atom_site_fract_y
    data_dict["_atom_site_fract_z"] = atom_site_fract_z
    data_dict["_atom_site_occupancy"] = atom_site_occupancy

    cif_str = str(cif_block)

    return cif_str


########################################


def read_cif(
    filename: str,
    workers: int = 1,
    keep_rare_gases: bool = False,
    keep_rare_earths: bool = False
) -> Tuple[List[Structure], int, int]:
    """
    Reads a cif file containing concatenated structures data and decode them using multiprocess.

    Parameters:
        filename (str):         Name of the input CIF file.
        
        workers (int):          Number of processes to use in parallel.
        
        keep_rare_gases (bool): Whether structures containing rare gases should be kept. 
                                Defaults to false.
        
        keep_rare_earths (bool): Whether structures containing f-block elements should be kept.
                                Defaults to False.
    
    Returns:
        List[Structure]: Decoded Structure objects in a list.
        Int: Number of structures containing rare gases discarded.
        Int: Number of structures containing rare earth elements discarded.
    """
    check_file_or_dir(filename, "file", allowed_formats="cif")
    check_type(workers, "workers", (int,))
    check_num_value(workers, "workers", ">", 0)
    check_type(keep_rare_gases, "keep_rare_gases", (bool,))
    check_type(keep_rare_earths, "keep_rare_earths", (bool,))

    struct_strings = extract_cif_from_file(filename)
    nbr_rare_gas_structs, nbr_rare_earth_structs = 0, 0

    if not keep_rare_gases:
        struct_strings, nbr_rare_gas_structs = discard_rare_gas_structures(struct_strings)
        print(f"{nbr_rare_gas_structs} rare gas structures ignored")

    if not keep_rare_earths:
        struct_strings, nbr_rare_earth_structs = discard_rare_earth_structures(struct_strings)
        print(f"{nbr_rare_earth_structs} rare earths structures ignored")

    nbr_struct = len(struct_strings)
    assert nbr_struct > 0, (
    "No structure data found in provided file (maybe they all have been discarded ?)"
    )
    chunksize = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    structs_list = list(process_map(
        cif_str_to_struct,
        struct_strings,
        max_workers=workers,
        chunksize=chunksize,
        desc="load and read data",
    ))

    return structs_list, nbr_rare_gas_structs, nbr_rare_earth_structs


########################################


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
    check_file_format(filename, allowed_formats="cif")
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))
    check_type(workers, "workers", (int,))
    check_num_value(workers, "workers", ">", 0)

    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    #def feed_args(structures, symprec, angle_tolerance) -> List[Tuple]:
    #    return [(struct, symprec, angle_tolerance) for struct in structures]

    encoded_cif = list(
        process_map(
            struct_to_cif_str,
            structures,
            max_workers=workers,
            chunksize=chunksize,
            desc="Writing data in cif format"
        )
    )

    with open(filename, "wt", encoding="utf-8") as out_file:
        out_file.write("\n".join(encoded_cif))


########################################
