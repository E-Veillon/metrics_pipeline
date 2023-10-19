#!/usr/bin/python
"""
A script that takes a .cif file containing crystal structures and returns a new .cif file with spacegroup symmetry calculated.
"""

##################################################
# SYSTEM I/O MODULES

import sys
from time import perf_counter
from typing import NamedTuple
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import re

##################################################
# OPTIMIZATION MODULES

# from functools import partial
# from itertools import zip_longest, count
from tqdm.contrib.concurrent import process_map

##################################################
# PYTHON MATERIALS GENOMICS MODULE

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

##################################################


def assert_args(args: NamedTuple):
    """
    Input arguments verification

    Args:
        args (NamedTuple): namespace of the parsed arguments.
    """

    assert args.filename.endswith(".cif") and args.output.endswith(
        ".cif"
    ), "some arguments formats are not supported, please only use .cif format"

    assert args.precision >= 0.0, "precision cannot be negative"

    assert (
        args.precision < 0.5
    ), "precision should stay lower than 0.5 to retain some reliability"

    assert args.angleprec >= 0.0, "angle tolerance cannot be negative"

    assert (
        args.angleprec < 20.0
    ), "angle tolerance should stay lower than 20.0 to retain some reliability"

    assert args.workers >= 1, "the number of workers cannot be negative or zero"


def has_rare_gas(formula: str) -> bool:
    """
    searches for rare gas symbols in structural formula

    Args:
        formula (NamedTuple): namespace of the parsed arguments.
    Returns
        bool: True if the formula contains rare gases, False otherwise.
    """

    return re.search(r"(He|Ne|Ar|Kr|Xe|Rn)", formula) is not None


def extract_from_cif_file(filename: str, keep_rare_gases: bool = False):
    """separates concatenated structures in a .cif file"""

    data_string_list = []
    full_structs_list = []
    containing_rare_gas = False
    rare_gas_structures = 0

    with open(filename, "r") as input_file:
        data = input_file.read()

    matcher = re.compile(r"^data.*?$(?=\ndata|\Z)", re.MULTILINE | re.DOTALL)
    structs = matcher.findall(data)

    rare_gas_structures = filter(has_rare_gas, structs)
    full_structs_list = filter(lambda x: not has_rare_gas(x), structs)

    return full_structs_list, rare_gas_structures


def cif_file_to_struct(filename: str) -> Structure:
    """Parses data from a .cif file and converts it into a pymatgen's Structure object (not recommended for files containing several structures)"""

    parsed_cif = CifParser(filename=filename)
    crystal_struct_list = parsed_cif.get_structures()
    crystal_struct = crystal_struct_list[0]
    symmetrized_struct = SpacegroupAnalyzer(
        crystal_struct, args.precision, args.angleprec
    ).get_symmetrized_structure()

    return symmetrized_struct


def cif_str_to_struct(string: str) -> Structure:
    """Parses data from a CIF formatted string and converts it into a pymatgen's Structure object"""

    parsed_str = CifParser.from_str(cif_string=string)
    crystal_struct_list = parsed_str.get_structures()
    crystal_struct = crystal_struct_list[0]
    symmetrized_struct = SpacegroupAnalyzer(
        crystal_struct, args.precision, args.angleprec
    ).get_symmetrized_structure()

    return symmetrized_struct


def structure_sort(groups_list: "list[list[Structure]]", keep_eliminated: bool = False):
    """Takes a list of structures grouped by equivalency and returns a sequence with unique structures only"""
    structures_kept = []
    structures_eliminated = []
    nbr_equivalent_structs = 0

    for group in groups_list:
        structures_kept.append(group[0])
        nbr_equivalent_structs += len(group[1:])

        if keep_eliminated:
            structures_eliminated.extend(group[1:])

    return structures_kept, structures_eliminated, nbr_equivalent_structs


def write_cif_string(struct: Structure):
    """Calculates Hermann-Mauguin's spacegroup in a pymatgen's Structure object and converts it into CIF formatted data string."""
    return str(
        CifWriter(
            struct=struct,
            symprec=args.precision,
            significant_figures=args.figures,
            angle_tolerance=args.angleprec,
        )
    )


def time_conversion(start: float, stop: float):
    """Converts time expressed in seconds to hours, minutes, seconds, milliseconds format."""

    total_seconds = stop - start
    seconds_left = int(total_seconds % 60)
    total_minutes = int(total_seconds // 60)
    minutes_left = int(total_minutes % 60)
    total_hours = int(total_minutes // 60)
    milliseconds = int((total_seconds % 60 - seconds_left) * 1000)

    return total_hours, minutes_left, seconds_left, milliseconds


if __name__ == "__main__":
    start = perf_counter()  # Initialisation de l'horloge

    # Gestion des arguments en ligne de commande

    prog_name = "Symmetrize_CIF"
    prog_description = "takes a .cif file containing crystal structures and returns a new .cif file with spacegroup symmetry calculated."
    help_format = ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, description=prog_description, formatter_class=help_format
    )

    parser.add_argument(
        "filename",
        type=str,
        help="name of input file to analyze",
        metavar="filename.cif",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="[filename]_out.cif",
        help="name of output file to return",
        metavar="output.cif",
    )
    parser.add_argument(
        "-p",
        "--precision",
        type=float,
        default=0.01,
        help="Fractional coordinates tolerance for symmetry finding",
        metavar="float",
    )
    parser.add_argument(
        "-a",
        "--angleprec",
        type=float,
        default=5.0,
        help="Angle tolerance for symmetry finding in degrees",
        metavar="float",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to create",
        metavar="int",
    )
    parser.add_argument(
        "-f",
        "--figures",
        type=int,
        default=8,
        help="Number of significant figures for written lattice parameters and sites coordinates",
        metavar="int",
    )
    parser.add_argument(
        "--keep_rare_gases",
        action="store_true",
        help="A flag to pass if data containing rare gases should not be automatically discarded",
    )
    parser.add_argument(
        "--keep_equivalent",
        action="store_true",
        help="A flag to pass if equivalent structures should be stored in a file instead of being discarded",
    )

    args = parser.parse_args()

    print(args)

    assert_args(args)

    if args.output == "[filename]_out.cif":
        args.output = args.filename.replace(".cif", "_out.cif")

    # Extraction des structures sous forme de strings

    struct_data, rare_gas_structures = extract_from_cif_file(
        filename=args.filename, keep_rare_gases=args.keep_rare_gases
    )

    if struct_data == []:
        sys.exit("Symmetrize_CIF: error: No structure data found in provided file")

    if not args.keep_rare_gases:
        print(f"{rare_gas_structures} structures containing rare gases discarded")

    # Obtention des structures PyMatGen à partir des données et calcul de la symétrie

    if len(struct_data) >= 200:
        chunksize = min(len(struct_data) // 100, 10)
    else:
        chunksize = 1

    print("Calculating spacegroup symmetry...")

    symmetrized_structs = process_map(
        cif_str_to_struct, struct_data, max_workers=args.workers, chunksize=chunksize
    )

    if symmetrized_structs == []:
        sys.exit("Symmetrize_CIF: error: No structure could be parsed from given data")

    # Comparaison des structures pour éliminer les doublons

    print("Comparing structures to find eventual equivalences...")

    matcher = StructureMatcher()
    grouped_structs = matcher.group_structures(s_list=symmetrized_structs)
    kept_structs, eliminated_structs, nbr_equiv = structure_sort(
        grouped_structs, keep_eliminated=args.keep_equivalent
    )

    print(f"{nbr_equiv} structures equivalent to one already calculated detected")

    if args.keep_equivalent:
        print(
            "The --keep_equivalent flag was passed, they will be stored in a separate file"
        )

    else:
        print("The --keep_equivalent flag wasn't passed, they will be discarded")

    #  Recalcul des symétries avec PyMatGen (pour prise en compte par CifWriter) et écriture du fichier de sortie

    print("Writing structures data into CIF format...")

    kept_cif_strings = process_map(
        write_cif_string, kept_structs, max_workers=args.workers, chunksize=chunksize
    )
    full_string = "\n".join(kept_cif_strings)

    with open(args.output, "wt") as out_file:
        out_file.write(full_string)

    if args.keep_equivalent:
        elim_cif_strings = process_map(
            write_cif_string,
            eliminated_structs,
            max_workers=args.workers,
            chunksize=chunksize,
        )
        full_string = "\n".join(elim_cif_strings)
        elim_structs_file = args.filename.replace(".cif", "_equiv.cif")

        with open(elim_structs_file, "wt") as elim_file:
            elim_file.write(full_string)

    # Calcul du temps total pris par la procédure

    stop = perf_counter()
    hours, minutes, seconds, milliseconds = time_conversion(start, stop)
    nbr_unique_structs = len(kept_structs)
    nbr_total_structs = rare_gas_structures + nbr_unique_structs + nbr_equiv
    print(" ")
    print("------------------------------")
    print(" ")
    print("SUMMARY OF THE CALCULATION")
    print(" ")
    print(f"{nbr_total_structs} structures detected in total, including:")
    print(f"- {rare_gas_structures} structures containing rare gases")
    print(f"- {nbr_equiv} structures equivalent to another one")
    print(f"- {nbr_unique_structs} unique structures")
    print(" ")
    print(f"Total calculation time: {hours}h {minutes}min {seconds}s {milliseconds}ms")
