#!/usr/bin/python
"""
A script to preprocess structures in a CIF file using pymatgen.
It filters out too bad and / or duplicated data, finds symmetry space group
and format the valid data in a way that will be more readable for further calculations.
"""

import os
import typing as typ
from datetime import datetime
import argparse as argp

# LOCAL IMPORTS
from . import _parse_input_args
from .utils import (
    check_type, check_num_value,
    check_file_format, check_file_or_dir,
    read_cif, write_cif, check_interatomic_distances,
    batch_symmetrizer, remove_equivalent, symmetrize_and_write_cif
)


########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_file",
        help="Path to the CIF file containing structure data to process.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Path to wanted CIF output file where processed structures "
            "will be stored. If not given, output file is written in input file "
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
        "-sk", "--special-keys",
        nargs="*",
        help=(
            "CIF Labels to store into structure properties, e.g. can be used "
            "to save structure identifiers attached to it throughout its manipulation "
            "as a python object."
        )
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
        help=(
            "Number of parallel processes to spawn for parallelized steps. "
            "If not given, If not given, default value is the 'max_workers' "
            "default value from tqdm.contrib.concurrent.process_map function."
        ),
        metavar="int",
    )
    parser.add_argument(
        "-S", "--sequential",
        action="store_true",
        help=(
            "Pass this flag to deactivate multiprocessing handling and switch to "
            "sequential computation. Generally, multiprocess is faster, but sometimes "
            "(e.g. when data is very big) it can softlock into resource distribution. "
            "Switch to more stable sequential computing if such problem  were to arise."
        )
    )
    args: argp.Namespace = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, typ.Any]) -> dict[str, typ.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    default_output = str(args_dict.get("input_file")).replace(".cif", "_out.cif")
    args_dict.setdefault("output", default_output)
    args_dict["output"] = os.path.realpath(args_dict["output"])
    args_dict.setdefault("no_rare_gas_check", False)
    args_dict.setdefault("no_rare_earth_check", False)
    args_dict.setdefault("no_dist_check", False)
    args_dict.setdefault("dist_tolerance", 0.5)
    args_dict.setdefault("no_symmetrization", False)
    args_dict.setdefault("symprec", 0.01)
    args_dict.setdefault("angleprec", 5.0)
    args_dict.setdefault("no_equiv_match", False)
    args_dict.setdefault("test_min_vol", False)
    args_dict.setdefault("sequential", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats="cif")
    check_file_format(args_dict.get("output"), allowed_formats="cif")
    check_type(args_dict.get("dist_tolerance"), "dist_tolerance", (float,))
    check_num_value(args_dict.get("dist_tolerance"), "dist_tolerance", ">", 0.0)
    check_type(args_dict.get("symprec"), "symprec", (float,))
    check_num_value(args_dict.get("symprec"), "symprec", ">=", 0.0)
    check_num_value(args_dict.get("symprec"), "symprec", "<=", 1.0)
    check_type(args_dict.get("angleprec"), "angleprec", (float,))
    check_num_value(args_dict.get("angleprec"), "angleprec", ">=", 0.0)
    check_num_value(args_dict.get("angleprec"), "angleprec", "<=", 90.0)

    args_dict["input_file"] = os.path.abspath(os.path.realpath(args_dict["input_file"]))
    args_dict["output"] = os.path.abspath(os.path.realpath(args_dict["output"]))

    if args_dict.get("special_keys") is not None:
        check_type(args_dict.get("special_keys"), "special_keys", (list,))
        for idx, elt in enumerate(args_dict["special_keys"]):
            check_type(elt, f"special_keys[{idx}]", (str,))

    if args_dict.get("workers") is not None:
        check_type(args_dict.get("workers"), "workers", (int,))
        check_num_value(args_dict.get("workers"), "workers", ">", 0)

    # Print final configuration
    print(" ")
    print("------------------------------")
    print(" ")
    print(" - I/O ARGUMENTS - ")
    print(f"INPUT FILE: {args_dict.get('input_file')}")
    print(f"OUTPUT FILE: {args_dict.get('output')}")
    print(f"SPECIAL KEYS: {'None' if args_dict.get('special_keys') is None else ''}")
    if args_dict.get("special_keys") is not None:
        keys_lst = '\n'.join(args_dict["special_keys"])
        print(f"{keys_lst}")
    print(" ")
    print("------------------------------")
    print(" ")
    print(" - ACTIVATED FEATURES - ")
    print(f"CHECK RARE GASES: {not args_dict.get('no_rare_gas_check')}")
    print(f"CHECK RARE EARTHS: {not args_dict.get('no_rare_earth_check')}")
    print(f"CHECK INTERATOMIC DISTANCES: {not args_dict.get('no_dist_check')}")
    print(
        f"* Distance tolerance: {args_dict.get('dist_tolerance')} Angstroms "
        f"{'(ignored)' if args_dict.get('no_dist_check') else ''}"
    )
    print(f"SYMMETRIZATION: {not args_dict.get('no_symmetrization')}")
    print(
        f"* Fractional coordinates tolerance: {args_dict.get('symprec')} "
        f"{'(ignored)' if args_dict.get('no_symmetrization') else ''}"
    )
    print(
        f"* Angles tolerance: {args_dict.get('angleprec')} degrees "
        f"{'(ignored)' if args_dict.get('no_symmetrization') else ''}"
    )
    print(f"STRUCTURE MATCHING: {not args_dict.get('no_equiv_match')}")
    print(f"IS SEQUENTIAL: {args_dict.get('sequential')}")
    print(
        "NUMBER OF WORKERS: "
        f"{'auto' if args_dict.get('workers') is None else args_dict.get('workers')} "
        f"{'(ignored)' if args_dict.get('sequential') else ''}"
    )
    print(" ")
    print("------------------------------")
    print(" ")

    return args_dict


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    A script to preprocess structures in a CIF file using pymatgen.
    It filters out too bad and / or duplicated data, finds symmetry space group
    and format the valid data in a way that will be more readable for further calculations.

    Args:
        standalone (bool):          Whether parsed script is used directly through
                                    command-line (stand-alone script) or in an external
                                    pipeline script.

        input_file (str|Path):      Path to the CIF file containing structure data to process.

        output (str|Path):          Path to wanted CIF output file where processed structures
                                    will be stored. If not given, output file is written in input
                                    file directory with the same name with a '_out' suffix added.

        no_rare_gas_check (bool):   Whether to disable elimination of structures containing rare
                                    gas elements. Defaults to False.

        no_rare_earth_check (bool): Whether to disable elimination of structures containing f-block
                                    elements. Defaults to False.

        no_dist_check (bool):       Whether to disable structures interatomic distances checking.

        dist_tolerance (float):     Tolerance for checking interatomic distances in Angstroms.
                                    Structures containing atoms that are closer than this value
                                    will be discarded from output. Defaults to 0.5 Angstroms.

        no_symmetrization (bool):   Whether to disable search of structures symmetry space groups.
                                    Defaults to False.

        symprec (float):            Fractional coordinates tolerance for symmetry finding.
                                    Defaults to 0.01.

        angleprec (float):          Angle tolerance for symmetry finding in degrees.
                                    Defaults to 5.0 degrees.

        no_equiv_match (bool):      Whether to disable structure matching and elimination of
                                    duplicates. Defaults to False.

        special_keys ([str]):       CIF Labels to store into structure properties, e.g.
                                    can be used to save structure identifiers attached to it
                                    throughout its manipulation as a python object.

        test_min_vol (bool):        A debug arg to assume unicity of unlikely structures having
                                    a volume under 1 Angström^3 without passing them into structure
                                    matching, which could cause the program to be softlocked. Only
                                    pass it if such problems were to arise when no_dist_check is
                                    set to True and no_equiv_match is set to False.

        workers (int):              Number of parallel processes to spawn for parallelized steps.
                                    If not given, If not given, default value is the 'max_workers'
                                    default value from tqdm.contrib.concurrent.process_map function.

        sequential (bool):          Pass this flag to deactivate multiprocessing handling and
                                    switch to sequential computation. Generally, multiprocess
                                    is faster, but sometimes (e.g. when data is very big) it
                                    can softlock into resource distribution. Switch to more
                                    stable sequential computing if such problem  were to arise.
    """
    start = datetime.now()
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    # Extraction des données CIF et conversion en structures
    structures, nbr_rare_gas_structs, nbr_rare_earth_structs = read_cif(
        filename=args["input_file"],
        keep_rare_gases=args["no_rare_gas_check"],
        keep_rare_earths=args["no_rare_earth_check"],
        special_keys=args.get("special_keys"),
        workers=args.get("workers"),
        sequential=args["sequential"]
    )

    nbr_loaded_structs = len(structures)
    nbr_total_structs  = nbr_loaded_structs + nbr_rare_gas_structs + nbr_rare_earth_structs
    assert nbr_loaded_structs > 0, "No structure could be parsed from given data"

    print(f"{nbr_total_structs} structures detected in total")

    if not args["no_rare_gas_check"]:
        print(f"{nbr_rare_gas_structs} structures containing rare gases were ignored")

    if not args["no_rare_earth_check"]:
        print(f"{nbr_rare_earth_structs} structures containing rare earths were ignored")

    print(f"{nbr_loaded_structs} structures are kept for further processing")

    # Vérification des distances interatomiques

    if not args["no_dist_check"]:
        structures, nbr_not_valid = check_interatomic_distances(
            structures, valid_tol=args["dist_tolerance"]
        )
        print(f"{nbr_not_valid} structures having too close atoms were discarded")

    # Calcul de la symétrie d'espace des structures

    """if args["no_symmetrization"]:
        symmetrized_structs = structures
    else:
        symmetrized_structs = list(filter(
            None,
            batch_symmetrizer(
                structures=structures,
                symprec=args["symprec"],
                angle_tolerance=args["angleprec"],
                workers=args["workers"],
                sequential=args["sequential"]
            )
        ))

        print(f"{len(symmetrized_structs)} structures were symmetrized")"""

    # Comparaison des structures pour éliminer les doublons

    kept_structs, nbr_equivalent, nbr_unmatched = remove_equivalent(
            structures=structures,
            workers=args.get("workers"),
            test_volume=args["test_min_vol"],
            keep_equivalent=args["no_equiv_match"],
            sequential=args["sequential"]
    )

    nbr_unique_structs = len(kept_structs)

    if not args["no_equiv_match"]:
        print(f"{nbr_unique_structs} unique structures detected")
        print(f"{nbr_equivalent} duplicates were discarded")

    if args["test_min_vol"]:
        print("--test-min-vol debug flag was passed:")
        print(f"{nbr_unmatched} structures were assumed unique.")

    # Symétrisation et écriture du fichier CIF épuré des structures indésirables
    if args["no_symmetrization"]:
        print("Symmetrization is deactivated, now writing CIF file...")
    else:
        print("Symmetrization is activated, now symmetrizing and writing CIF file...")

    symmetrize_and_write_cif(
        filename=args["output"],
        structures=kept_structs,
        symmetrize=(not args["no_symmetrization"]),
        symprec=args["symprec"],
        angleprec=args["angleprec"],
        special_keys=args.get("special_keys"),
        workers=args.get("workers"),
        sequential=args["sequential"]
    )

    # Calcul du temps total pris par la procédure

    stop = datetime.now()

    print("\n------------------------------")
    print("\nSUMMARY OF THE CALCULATION")
    print(f"{nbr_total_structs} structures detected in total, including:")
    print(f"- {nbr_unique_structs} unique structure(s)")

    if not args["no_rare_gas_check"]:
        print(f"- {nbr_rare_gas_structs} structure(s) containing rare gases")

    if not args["no_rare_earth_check"]:
        print(f"- {nbr_rare_earth_structs} structure(s) containing rare earths")

    if not args["no_dist_check"]:
        print(f"- {nbr_not_valid} structure(s) with too small interatomic distances")

    if not args["no_equiv_match"]:
        print(f"- {nbr_equivalent} structure(s) that are duplicates")

    if args["test_min_vol"]:
        print(
            f"- {nbr_unmatched} structure(s) assumed unique "
            f"because of {'its' if nbr_unmatched < 2 else 'their'} unphysical volume"
        )
    print(f"\nOutput results written in {args['output']}")
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
