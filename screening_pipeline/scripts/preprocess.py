#!/usr/bin/python
"""
A script to preprocess structures in a CIF file using pymatgen.
It filters out too bad and / or duplicated data, finds symmetry space group
and format the valid data in a way that will be more readable for further calculations.
"""

from datetime import datetime
import argparse

# LOCAL IMPORTS
from screening_pipeline.utils import (
    check_file_format, check_file_or_dir,
    read_cif, write_cif, check_interatomic_distances,
    batch_symmetrizer, remove_equivalent
)


########################################


def _assert_args(args: argparse.Namespace):
    """Asserting input arguments validity."""

    print(" - I/O ARGUMENTS - ")
    print(f"INPUT FILE: {args.input_file}")
    print(f"OUTPUT FILE: {args.output}")

    check_file_or_dir(args.input_file, "file", allowed_formats="cif")
    check_file_format(args.output, allowed_formats="cif")

    print(" ")
    print("------------------------------")
    print(" ")
    print(" - ACTIVATED FEATURES - ")
    print(f"CHECK RARE GASES: {not args.no_rare_gas_check}")
    print(f"CHECK RARE EARTHS: {not args.no_rare_earth_check}")
    print(f"CHECK INTERATOMIC DISTANCES: {not args.no_dist_check}")
    print(
        f"* Distance tolerance: {args.dist_tolerance} Angstroms "
        f"{'(ignored)' if args.no_dist_check else ''}"
    )

    if not args.no_dist_check:
        assert args.dist_tolerance > 0.0, (
            "Interatomic distance tolerance must be strictly positive."
        )
    print(f"SYMMETRIZATION: {not args.no_symmetrization}")
    print(
        f"* Fractional coordinates tolerance: {args.symprec} "
        f"{'(ignored)' if args.no_symmetrization else ''}"
    )

    if not args.no_symmetrization:
        assert (0.0 <= args.symprec <= 0.5), (
        "Fractional coordinates tolerance must be between 0.0 and 0.5 "
        "to retain some reliability."
        )

    print(
        f"* Angles tolerance: {args.angleprec} degrees "
        f"{'(ignored)' if args.no_symmetrization else ''}"
    )

    if not args.no_symmetrization:
        assert (0.0 <= args.angleprec <= 30.0), (
            "Angles tolerance must be between 0.0 and 20.0 degrees to retain some reliability."
        )
    print(f"STRUCTURE MATCHING: {not args.no_equiv_match}")
    print(f"NUMBER OF WORKERS: {args.workers}")

    assert args.workers >= 1, "The number of workers must be positive."

    print(" ")
    print("------------------------------")
    print(" ")


########################################


def main() -> None:
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "preprocess.py"
    prog_description = (
        "A script to preprocess structures in a CIF file using pymatgen. "
        "It filters out too bad and duplicated data, finds symmetry space group "
        "and format the valid data in a way that will be more readable "
        "for further calculations."
    )

    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=prog_description
    )

    parser.add_argument(
        "input_file",
        type=str,
        help="Name or path to the CIF file containing structure data to process.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help=(
            "Name or path to wanted CIF output file where processed structures "
            "will be stored. If not specified, output file is written in input file "
            "directory with the same name with a '_out' suffix added."
        ),
    )
    parser.add_argument(
        "--no-rare-gas-check",
        action="store_true",
        help="A flag to disable elimination of structures containing rare gas elements.",
        dest="no_rare_gas_check"
    )
    parser.add_argument(
        "--no-rare-earth-check",
        action="store_true",
        help="A flag to disable elimination of structures containing f-block elements.",
        dest="no_rare_earth_check"
    )
    parser.add_argument(
        "--no-dist-check",
        action="store_true",
        help="A flag to disable structures interatomic distances checking.",
        dest="no_dist_check"
    )
    parser.add_argument(
        "-d",
        "--dist-tolerance",
        type=float,
        default=0.5,
        help=(
            "Tolerance for checking interatomic distances in Angstroms. "
            "Structures containing atoms that are closer than this value will be discarded "
            "(Default: %(default)s Angstroms)."
        ),
        metavar="float",
        dest="dist_tolerance"
    )
    parser.add_argument(
        "--no-symmetrization",
        action="store_true",
        help="A flag to disable search of structures symmetry space groups.",
        dest="no_symmetrization"
    )
    parser.add_argument(
        "-s",
        "--symprec",
        type=float,
        default=0.01,
        help="Fractional coordinates tolerance for symmetry finding (Default: %(default)s).",
        metavar="float",
    )
    parser.add_argument(
        "-a",
        "--angleprec",
        type=float,
        default=5.0,
        help="Angle tolerance for symmetry finding in degrees (Default: %(default)s degrees).",
        metavar="float",
    )
    parser.add_argument(
        "--no-equiv-match",
        action="store_true",
        help="A flag to disable structure matching and elimination of duplicates.",
        dest="no_equiv_match"
    )
    parser.add_argument(
        "--test-min-vol",
        action="store_true",
        help=(
            "A debug flag to assume unicity of unlikely structures having a volume under "
            "1 Angström^3 without passing them into structure matching, which could cause "
            "the program to be softlocked. Only pass it if such problems were to arise when "
            "--no-dist-check flag is passed and not --no-equiv-match flag."
        )
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to spawn (Default: %(default)s).",
        metavar="int",
    )

    args: argparse.Namespace = parser.parse_args()

    if args.output is None:
        args.output = args.input_file.replace(".cif", "_out.cif")

    _assert_args(args)


    # MAIN BLOCK

    # Extraction des données CIF et conversion en structures

    structures, nbr_rare_gas_structs, nbr_rare_earth_structs = read_cif(
        filename=args.input_file,
        workers=args.workers,
        keep_rare_gases=args.no_rare_gas_check,
        keep_rare_earths=args.no_rare_earth_check
    )

    nbr_loaded_structs = len(structures)
    nbr_total_structs  = nbr_loaded_structs + nbr_rare_gas_structs + nbr_rare_earth_structs
    assert nbr_loaded_structs > 0, "No structure could be parsed from given data"

    print(f"{nbr_total_structs} structures detected in total")

    if not args.no_rare_gas_check:
        print(f"{nbr_rare_gas_structs} structures containing rare gases were ignored")

    if not args.no_rare_earth_check:
        print(f"{nbr_rare_earth_structs} structures containing rare earths were ignored")

    print(f"{nbr_loaded_structs} structures are kept for further processing")

    # Vérification des distances interatomiques

    if not args.no_dist_check:
        structures, nbr_not_valid = check_interatomic_distances(
            structures, valid_tol=args.dist_tolerance
        )
        print(f"{nbr_not_valid} structures having too close atoms were discarded")

    # Calcul de la symétrie d'espace des structures

    if args.no_symmetrization:
        symmetrized_structs = structures
    else:
        symmetrized_structs = list(filter(
            None,
            batch_symmetrizer(
                structures=structures,
                symprec=args.symprec,
                angle_tolerance=args.angleprec,
                workers=args.workers
            )
        ))

        print(f"{len(symmetrized_structs)} structures were symmetrized")

    # Comparaison des structures pour éliminer les doublons

    kept_structs, nbr_equivalent, nbr_unmatched = remove_equivalent(
            structures=symmetrized_structs,
            workers=args.workers,
            test_volume=args.test_min_vol,
            keep_equivalent=args.no_equiv_match
    )

    nbr_unique_structs = len(kept_structs)

    if not args.no_equiv_match:
        print(f"{nbr_unique_structs} unique structures detected")
        print(f"{nbr_equivalent} duplicates were discarded")

    if args.test_min_vol:
        print("--test-min-vol debug flag was passed:")
        print(f"{nbr_unmatched} structures were assumed unique.")

    # Ecriture du fichier CIF symétrisé et épuré des structures indésirables
    write_cif(
        filename=args.output,
        structures=kept_structs,
        workers=args.workers,
    )

    # Calcul du temps total pris par la procédure

    stop = datetime.now()

    print(" ")
    print("------------------------------")
    print(" ")
    print("SUMMARY OF THE CALCULATION")
    print(" ")
    print(f"{nbr_total_structs} structures detected in total, including:")
    print(f"- {nbr_unique_structs} unique structure(s)")

    if not args.no_rare_gas_check:
        print(f"- {nbr_rare_gas_structs} structure(s) containing rare gases")

    if not args.no_rare_earth_check:
        print(f"- {nbr_rare_earth_structs} structure(s) containing rare earths")

    if not args.no_dist_check:
        print(f"- {nbr_not_valid} structure(s) with too small interatomic distances")

    if not args.no_equiv_match:
        print(f"- {nbr_equivalent} structure(s) that are duplicates")

    if args.test_min_vol:
        print(
            f"- {nbr_unmatched} structure(s) assumed unique "
            f"because of {'its' if nbr_unmatched < 2 else 'their'} unphysical volume"
        )

    print(" ")
    print(f"Output results written in '{args.output}'")
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
