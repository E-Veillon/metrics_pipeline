#!/usr/bin/python
"""
A script that takes a .cif file containing crystal structures and returns a new .cif file with spacegroup symmetry calculated.
"""

##################################################
# SYSTEM I/O MODULES

from datetime import datetime
from typing import NamedTuple
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def assert_args(args: NamedTuple):
    """
    Input arguments verification

    Args:
        args (NamedTuple): namespace of the parsed arguments.
    """

    assert (
        args.filename.endswith(".cif")
        and args.output.endswith(".cif")
        and (args.equivalent is None or args.equivalent.endswith(".cif"))
    ), "some arguments formats are not supported, please only use .cif format"

    assert (
        0.0 <= args.precision <= 0.5
    ), "precision must be between 0.0 and 0.5 Angstrom"

    assert (
        0.0 <= args.angleprec <= 20.0
    ), "angle tolerance must be between 0.0 and 20.0 degree"

    assert args.workers >= 1, "the number of workers cannot be negative or zero"


def main():
    start = datetime.now()

    # Gestion des arguments en ligne de commande

    prog_name = "symmetrize"
    prog_description = "takes a .cif file containing crystal structures and returns a new .cif file with spacegroup symmetry calculated."
    help_format = ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, description=prog_description, formatter_class=help_format
    )

    parser.add_argument(
        "filename",
        type=str,
        help="name of the file containing the input structures",
        default="filename.cif",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="[filename]_out.cif",
        help="name of output file containing the unique structures",
    )
    parser.add_argument(
        "-e",
        "--equivalent",
        type=str,
        default=None,
        help="name of the output file containing the duplicated structures",
    )
    parser.add_argument(
        "-p",
        "--precision",
        type=float,
        default=0.5,
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
        "--keep-rare-gases",
        action="store_true",
        help="A flag to pass if data containing rare gases should not be automatically discarded",
    )
    parser.add_argument(
        "--keep-equivalent",
        action="store_true",
        help="A flag to pass if data containing rare gases should not be automatically discarded",
    )

    args = parser.parse_args()

    print(args)

    assert_args(args)

    from screening_pipeline.utils import (
        has_rare_gas,
        read_cif,
        write_cif,
        remove_equivalent,
    )

    if args.output == "[filename]_out.cif":
        args.output = args.filename.replace(".cif", "_out.cif")

    # Extraction des structures sous forme de strings

    symmetrized_structs = read_cif(
        args.filename,
        symprec=args.precision,
        angle_tolerance=args.angleprec,
        workers=args.workers,
        keep_rare_gases=args.keep_rare_gases,
    )

    assert len(symmetrized_structs) > 0, "No structure could be parsed from given data"

    print(f"{len(symmetrized_structs)} structures loaded")

    # Comparaison des structures pour éliminer les doublons

    if args.keep_equivalent:
        kept_structs = symmetrized_structs
        duplicated_struct = []
    else:
        kept_structs, duplicated_struct = remove_equivalent(
            symmetrized_structs, workers=args.workers, keep_equivalent=False
        )

        print(f"{len(kept_structs)} unique structures detected")

    #  Recalcul des symétries avec PyMatGen (pour prise en compte par CifWriter) et écriture du fichier de sortie

    write_cif(
        args.output,
        kept_structs,
        symprec=args.precision,
        angle_tolerance=args.angleprec,
        workers=args.workers,
    )

    if args.equivalent is not None:
        write_cif(
            args.equivalent,
            duplicated_struct,
            symprec=args.precision,
            angle_tolerance=args.angleprec,
            workers=args.workers,
        )

    # Calcul du temps total pris par la procédure

    stop = datetime.now()
    count_unique = len(kept_structs)
    count_duplicated = len(duplicated_struct)
    print(" ")
    print("------------------------------")
    print(" ")
    print("SUMMARY OF THE CALCULATION")
    print(" ")
    print(f"{count_duplicated+count_unique} structures detected in total, including:")
    print(f"- {count_duplicated} duplicated structures")
    print(f"- {count_unique} unique structures")
    print(" ")
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
