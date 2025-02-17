#!/usr/bin/python
"""Process CSV data according to its 'cif' column."""

import os
import argparse as argp
from tqdm import tqdm

import pandas as pd
import numpy as np

from pymatgen.io.cif import CifParser, CifWriter

from .utils import (
    check_interatomic_distances,
    remove_equivalent,
    discard_rare_gas_structures,
    discard_rare_earth_structures,
    check_file_or_dir, check_file_format, check_num_value
)


def _parse_input_args() -> argp.Namespace:
    parser = argp.ArgumentParser(description=__doc__)
    parser.add_argument("infile",help="CSV file containing a 'cif' data column to process.")
    parser.add_argument("-o","--output",help="Path to write processed CSV file.")
    parser.add_argument("-w","--workers",type=int,help="Number of parallel processes to spawn.")
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

    check_file_or_dir(args.infile, "file", allowed_formats="csv")
    
    if args.output is None:
        args.output = args.infile.replace(".csv","_parsed.csv")
    else:
        check_file_format(args.output, allowed_formats="csv")
    
    if args.workers is not None:
        check_num_value(args.workers, "args.workers", ">", 0)
    
    check_num_value(args.dist_tolerance, "--dist-tolerance", ">", 0)

    return args


def main(args: argp.Namespace|None = None, **kwargs) -> None:
    df = pd.read_csv(args.infile)
    ids: list[int] = df["material_id"].tolist()
    cifs: list[str] = df["cif"].tolist()
    step_sep = "\n------------------------------\n"

    # Rename CIFs headers with their 'material_id'
    new_cifs = []
    for id, cif in tqdm(list(zip(ids, cifs)),desc="Renaming CIFs"):
        split_cif = cif.splitlines()
        split_cif[0] = "data_oqmd-" + str(id)
        new_cif = "\n".join(split_cif)
        new_cifs.append(new_cif)
        df.at[ids.index(id), "cif"] = new_cif
    new_data = new_cifs
    print(step_sep)

    # Process CIF strings
    if not args.no_rare_gas_check:
        print("Searching for rare gases...")
        no_rg_cifs, nbr_rg = discard_rare_gas_structures(new_cifs)
        new_data = no_rg_cifs
        print(f"Number of structures containing rare gases: {nbr_rg}")
        print(f"Structures left: {len(no_rg_cifs)}")
        print(step_sep)

    if not args.no_rare_earth_check:
        print("Searching for rare earths...")
        no_rg_re_cifs, nbr_re = discard_rare_earth_structures(new_data)
        new_data = no_rg_re_cifs
        print(f"Number of structures containing rare earths: {nbr_re}")
        print(f"Structures left: {len(no_rg_re_cifs)}")
        print(step_sep)

    # Structurize kept CIFs and save their id
    structs = []
    for cif in tqdm(new_data, desc="Structurize CIFs"):
        id = int(cif.splitlines()[0][10:])
        struct = CifParser.from_str(cif).parse_structures(primitive=False, on_error='ignore')[0]
        struct.properties["id"] = id
        structs.append(struct)
    new_data = structs
    print(step_sep)

    if not args.no_dist_check:
        # Process structures
        print(f"Searching for non-valid structures (atom pairs closer than {args.dist_tolerance}A)...")
        valid_structs, nbr_no_valid = check_interatomic_distances(new_data, valid_tol=args.dist_tolerance)
        new_data = valid_structs
        print(f"Number of not valid structures: {nbr_no_valid}")
        print(f"Structures left: {len(valid_structs)}")
        print(step_sep)

    if not args.no_equiv_match:
        print("Searching for duplicates...")
        val_uniq_structs, nbr_dupl, _ = remove_equivalent(new_data, workers=args.workers)
        new_data = val_uniq_structs
        print(f"Number of duplicate structures: {nbr_dupl}")
        print(f"Structures left: {len(val_uniq_structs)}")
        print(step_sep)

    print("Build new DataFrame with parsed data...")
    # Get 'material_id's of all kept data
    kept_ids = set([struct.properties["id"] for struct in new_data])
    print(f"Number of kept ids: {len(kept_ids)}")

    # Build a new DataFrame with valid structures data only from the old DataFrame
    def by_Series_idx(row_tuple: tuple[int, pd.Series]):
        return row_tuple[0]

    kept_rows = list(filter(lambda row: row[1].material_id in kept_ids, df.iterrows()))
    print(f"Number of kept data rows: {len(kept_rows)}")
    sorted_rows = [row[1] for row in sorted(kept_rows, key=by_Series_idx)]
    print(f"Number of kept data rows sorted: {len(sorted_rows)}")
    parsed_df = pd.DataFrame(data=sorted_rows)
    print(f"DataFrame built, {len(parsed_df)=}.")
    print(step_sep)

    # Write new DataFrame to CSV file
    print("Write new DataFrame to CSV file...")
    parsed_df.to_csv(args.output)
    print(f"Parsed data written in '{args.output}'.")


if __name__ == "__main__":
    args = _parse_input_args()
    main(args)
