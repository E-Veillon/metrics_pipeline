#!/usr/bin/python
"""Process CSV data according to its 'cif' column."""

import os
import typing as tp
import argparse as ap
from tqdm import tqdm

import pandas as pd

from src.utils import (
    parse_input_args, check_type, check_num_value,
    discard_rare_gas_structures, discard_rare_earth_structures
)
from src.genmat_io import CIFFile, check_file_or_dir, check_file_format
from src.metrics import StructValidity, Unicity


def _get_command_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument("infile",help="CSV file containing a 'cif' data column to process.")
    parser.add_argument("-o","--output",help="Path to write processed CSV file.")
    parser.add_argument(
        "--id-key", default="material_id",
        help=(
            "Name of the column containing structures unique identifier. "
            "Defaults to %(default)s."
        )
    )
    parser.add_argument(
        "-w","--workers",
        type=int,
        help=(
            "Number of parallel processes to spawn for parallelized steps. "
            "If not given, If not given, default value is the 'max_workers' "
            "default value from tqdm.contrib.concurrent.process_map function."
        )
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
    args = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    default_output = str(args_dict.get("infile")).replace(".csv","_parsed.csv")
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("id_key", "material_id")
    args_dict.setdefault("no_rare_gas_check", False)
    args_dict.setdefault("no_rare_earth_check", False)
    args_dict.setdefault("no_dist_check", False)
    args_dict.setdefault("dist_tolerance", 0.5)
    args_dict.setdefault("no_equiv_match", False)
    args_dict.setdefault("test_min_vol", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("infile"), "file", allowed_formats="csv")
    check_file_format(args_dict.get("output"), allowed_formats="csv")
    check_type(args_dict["id_key"], "id_key", (str,))

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">", 0)
    
    check_type(args_dict["dist_tolerance"], "dist_tolerance", (float,))
    check_num_value(args_dict["dist_tolerance"], "dist_tolerance", ">", 0.0)

    return args_dict


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    Process CSV data according to its 'cif' column.
    
    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    infile: str | Path
        CSV file containing a 'cif' data column to process.

    output: str | Path
        Path to write processed CSV file.

    id_key: str, optional
        Name of the column containing structures unique identifier.
        Defaults to %(default)s.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    no_rare_gas_check: bool
        Whether to disable elimination of structures containing rare gas elements.
        Defaults to False.

    no_rare_earth_check: bool
        Whether to disable elimination of structures containing f-block elements.
        Defaults to False.

    no_dist_check: bool
        Whether to disable structures interatomic distances checking.

    dist_tolerance: float
        Tolerance for checking interatomic distances in Angstroms. Structures containing
        atoms that are closer than this value will be discarded from output.
        Defaults to 0.5 Angstroms.

    no_equiv_match: bool
        Whether to disable structure matching and elimination of duplicates.
        Defaults to False.

    test_min_vol: bool
        A debug arg to assume unicity of unlikely structures having a volume under 1 Angström^3
        without passing them into structure matching, which could cause the program to be
        softlocked. Only pass it if such problems were to arise when no_dist_check is set to
        True and no_equiv_match is set to False.
    """
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    df = pd.read_csv(args["infile"])
    ids: list[int] = df[args["id_key"]].tolist()
    cifs: list[str] = df["cif"].tolist()
    step_sep = "\n------------------------------\n"

    # Rename CIFs headers with their 'id_key'
    new_cifs = []
    for id, cif in tqdm(list(zip(ids, cifs)),desc="Renaming CIFs"):
        split_cif = cif.split(sep="\n", maxsplit=1)
        split_cif[0] = "data_" + str(id)
        new_cif = "\n".join(split_cif)
        new_cifs.append(new_cif)
        df.at[ids.index(id), "cif"] = new_cif
    new_data = new_cifs
    print(step_sep)

    # Process CIF strings
    if not args.get("no_rare_gas_check"):
        print("Searching for rare gases...")
        no_rg_cifs, nbr_rg = discard_rare_gas_structures(new_data)
        new_data = no_rg_cifs
        print(f"Number of structures containing rare gases: {nbr_rg}")
        print(f"Structures left: {len(no_rg_cifs)}")
        print(step_sep)

    if not args.get("no_rare_earth_check"):
        print("Searching for rare earths...")
        no_rg_re_cifs, nbr_re = discard_rare_earth_structures(new_data)
        new_data = no_rg_re_cifs
        print(f"Number of structures containing rare earths: {nbr_re}")
        print(f"Structures left: {len(no_rg_re_cifs)}")
        print(step_sep)

    # Structurize kept CIFs and save their id
    structs, _ = CIFFile(
        new_data, special_keys=["header"], workers=args.get("workers")
    ).parse_structures()

    new_data = structs
    print(step_sep)

    if not args.get("no_dist_check"):
        # Process structures
        print(
            "Searching for non-valid structures "
            f"(atom pairs closer than {args.get('dist_tolerance')}A)..."
        )
        validity = StructValidity(new_data, valid_tol=args["dist_tolerance"])
        valid_structs, nbr_no_valid = validity.valid_structs, len(validity.invalid_structs)
        new_data = valid_structs
        print(f"Number of not valid structures: {nbr_no_valid}")
        print(f"Structures left: {len(valid_structs)}")
        print(step_sep)

    if not args.get("no_equiv_match"):
        print("Searching for duplicates...")
        unicity = Unicity(new_data, workers=args.get("workers"))
        val_uniq_structs, nbr_dupl = unicity.unique_structs, len(unicity.duplicate_structs)
        new_data = val_uniq_structs
        print(f"Number of duplicate structures: {nbr_dupl}")
        print(f"Structures left: {len(val_uniq_structs)}")
        print(step_sep)

    print("Build new DataFrame with parsed data...")
    # Get IDs of all kept data
    kept_ids = set([struct.properties["header"] for struct in new_data])
    print(f"Number of kept ids: {len(kept_ids)}")

    # Build a new DataFrame with valid structures data only from the old DataFrame
    def by_Series_idx(row_tuple: tuple[int, pd.Series]) -> int:
        return row_tuple[0]

    kept_rows = list(
        filter(lambda row: getattr(row[1], args["id_key"]) in kept_ids, df.iterrows())
    )
    print(f"Number of kept data rows: {len(kept_rows)}")
    sorted_rows = [row[1] for row in sorted(kept_rows, key=by_Series_idx)] # type: ignore
    print(f"Number of kept data rows sorted: {len(sorted_rows)}")
    parsed_df = pd.DataFrame(data=sorted_rows)
    print(f"DataFrame built, {len(parsed_df)=}.")
    print(step_sep)

    # Write new DataFrame to CSV file
    print("Write new DataFrame to CSV file...")
    parsed_df.to_csv(args.get("output"))
    print(f"Parsed data written in {args.get('output')}.")


if __name__ == "__main__":
    main()
