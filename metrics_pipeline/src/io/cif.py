"""I/O operations for concatenated CIF data files."""

import os
import re
import functools as ft

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser, CifWriter

from .io_base import PathLike, check_file_or_dir, check_file_format
from src.utils import VisualIterator

# TODO: Split responsibilities better, e.g. define compound convenience functions
# outside the package.
from src.utils.periodic_table import (
    discard_rare_gas_structures, discard_rare_earth_structures
)


def extract_cif_from_file(filename: PathLike) -> list[str]:
    """
    Reads and separate concatenated structures from a cif file.

    Parameters
    ----------
    filename: str
        Path to a cif file.

    Returns
    -------
    list[str]
        List of CIF strings, each one representing a single structure.
    """
    check_file_or_dir(filename, "file", allowed_formats="cif")

    # load file
    with open(filename, "rt", encoding="utf-8") as input_file:
        data = input_file.read()

    # match structures data with regular expression
    matcher = re.compile(r"^data.*?$(?=\ndata|\Z)", re.MULTILINE | re.DOTALL)
    cif_str_structs = [str(cif_struct).strip() for cif_struct in matcher.findall(data)]

    return cif_str_structs


def cif_str_to_struct(
    cif_str: str, special_keys: list[str] | None = None
) -> Structure:
    """
    Convert data from a cif formatted string to a pymatgen Structure object.

    Parameters
    ----------
    cif_str: str
        The string of a structure encoded in the cif format.

    special_keys: list[str], optional
        CIF Labels to store into structure properties, e.g. can be used to save
        structure identifiers attached to it throughout its manipulation as a
        python object. Pass "header" to save the original header.

    Returns
    -------
    Structure
        A Pymatgen Structure object.
    """
    parser = CifParser.from_str(cif_string=cif_str)
    structure = parser.parse_structures(primitive=False)[0]

    if special_keys is not None:
        props = dict.fromkeys(special_keys)
        cif_dict: dict = parser.as_dict().popitem()[1]
        for key in special_keys:
            if key == "header":
                props[key] = parser._cif.data.popitem()[1].header
            else:
                props[key] = cif_dict.get(key)

        structure.properties.update(props)

    return structure


def struct_to_sym_cif_str(
    structure: Structure,
    significant_figures: int = 8,
    symprec: float|None = 0.01,
    angleprec: float = 5.0,
    special_keys: list[str]|None = None
) -> str:
    """
    Symmetrize and converts a pymatgen Structure object into a CIF formatted string.

    Parameters:
        structure (Structure):      Structure object to convert.

        significant_figures (int):  Number of decimal places to keep for atomic positions.

        symprec (float|None):       Fractional position tolerance to find symmetry.
                                    If set to None, symmetry finding is disabled.
                                    Defaults to 0.01.

        angleprec (float):          Angle tolerance to find symmetry. Defaults to 5 degrees.

        special_keys ([str]):       CIF Labels to store into structure properties, e.g.
                                    can be used to save structure identifiers attached to it
                                    throughout its manipulation as a python object.

    Returns:
        str: CIF formatted string.
    """
    writer = CifWriter(
            struct=structure,
            symprec=symprec,
            significant_figures=significant_figures,
            angle_tolerance=angleprec
        )

    if special_keys is not None:
        final_cif = []
        auto_cif = str(writer).splitlines(keepends=True)

        # Remove comment lines
        auto_cif = list(filter(lambda line: not line.startswith("#"), auto_cif))

        # Recovery of the saved header
        if "header" in special_keys:
            final_cif.append(f"data_{structure.properties.pop('header')}\n")
        else:
            final_cif.append(auto_cif[0])

        # Recovery of saved attributes
        for key in special_keys:
            if key == "header":
                continue
            final_cif.append(f"{key}   {structure.properties.get(key)}\n")

        # Addition of the generated CIF data
        final_cif.extend(auto_cif[1:])
        final_cif = "".join(final_cif)

    else:
        final_cif = str(writer)

    return final_cif


def read_cif(
    filename: PathLike,
    keep_rare_gases: bool = False,
    keep_rare_earths: bool = False,
    special_keys: list[str] | None = None,
    workers: int | None = None,
) -> tuple[list[Structure], int, int]:
    """
    Read a cif file containing concatenated structures data and parse them using multiprocess.

    Parameters
    ----------
    filename: str | Path
        Path to the input CIF file.

    keep_rare_gases: bool
        Whether structures containing rare gases should be kept. Defaults to false.

    keep_rare_earths (bool
        Whether structures containing f-block elements should be kept. Defaults to False.

    special_keys: list[str], optional
        CIF Labels to store into structure properties, e.g. can be used to save structure
        identifiers attached to it throughout its manipulation as a python object.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    Returns
    -------
    List[Structure]
        Decoded Structure objects in a list.
    Int
        Number of structures containing rare gases discarded.
    Int
        Number of structures containing rare earth elements discarded.
    """
    check_file_or_dir(filename, "file", allowed_formats="cif")

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

    if special_keys is not None:
        partial_fn = ft.partial(cif_str_to_struct, special_keys=special_keys)
    else:
        partial_fn = cif_str_to_struct

    if workers is not None and workers == 0:
        data_list = [
            partial_fn(struct_string)
            for struct_string in VisualIterator(
                struct_strings, desc="Read and load data", unit="structures", percent=True
            )
        ]
    else:
        data_list = list(process_map(
            partial_fn,
            struct_strings,
            max_workers=workers,
            chunksize=chunksize,
            desc="load and read data",
        ))
    structs_list: list[Structure] = list(filter(lambda data: isinstance(data, Structure), data_list))
    unparsed_cifs: list[str] = list(filter(lambda data: isinstance(data, str), data_list))
    if unparsed_cifs:
        print(
            f"WARNING: {len(unparsed_cifs)} data could not be parsed into "
            "Structure objects, they will be removed from processing and "
            "written as-is in a 'pmg_unparsed.cif' file."
        )
        with open("pmg_unparsed.cif","wt") as fp:
            fp.write("\n".join(unparsed_cifs))

    return structs_list, nbr_rare_gas_structs, nbr_rare_earth_structs


def symmetrize_and_write_cif(
    filename: str,
    structures: list[Structure],
    significant_figures: int = 8,
    symmetrize: bool = True,
    symprec: float = 0.01,
    angleprec: float = 5.0,
    special_keys: list[str]|None = None,
    workers: int|None = None
) -> None:
    """
    Symmetrize and encode multiple structures in CIF format and write them in a file
    using multiprocess.

    Parameters
    ----------
    filename (str
        Name of the input file.

    structures (List[Structure]
        The structures to encode.

    significant_figures (int
        Number of decimal places to keep for atomic positions.

    symmetrize (bool
        Whether to search for structures spacegroup symmetry and refine their atomic positions
        according to found symmetry before encoding them. Defaults to True.

    symprec (float
        Fractional position tolerance to find symmetry. Defaults to 0.01.

    angleprec (float
        Angle tolerance to find symmetry. Defaults to 5 degrees.

    special_keys ([str]
        CIF Labels to store into structure properties, e.g. can be used to save structure
        identifiers attached to it throughout its manipulation as a python object.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.
    """
    check_file_format(filename, allowed_formats="cif")

    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    if symmetrize:
        partial_fn = ft.partial(
            struct_to_sym_cif_str,
            significant_figures=significant_figures,
            symprec=symprec,
            angleprec=angleprec,
            special_keys=special_keys
        )
    else:
        partial_fn = ft.partial(
            struct_to_sym_cif_str,
            significant_figures=significant_figures,
            symprec=None,
            special_keys=special_keys
        )

    if workers == 0:
        encoded_cif = [
            partial_fn(struct)
            for struct in VisualIterator(structures, desc="Writing data in cif format")
        ]

    else:
        encoded_cif = list(
            process_map(
                partial_fn,
                structures,
                max_workers=workers,
                chunksize=chunksize,
                desc="Writing data in cif format"
            )
        )
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "wt", encoding="utf-8") as out_file:
        out_file.write("\n".join(encoded_cif))