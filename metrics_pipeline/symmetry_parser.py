"""
Find symmetry of CIF structures, refine positions according to symmetry
and store them in different files depending on their symmetry class.
"""

import os.path as path
from enum import Enum
import argparse as argp
import typing as tp
import functools as ft

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from src.utils import parse_input_args, check_type, check_num_value, VisualIterator, PG_TO_SYSTEM
from src.io import CIFFile, check_file_or_dir


class SymmetryClass(Enum):
    FAMILY = "family"
    SYSTEM = "system"
    POINTGROUP = "pointgroup"
    SPACEGROUP = "spacegroup"


def _get_command_line_args() -> argp.Namespace:
    """Handle CLI arguments."""
    parser = argp.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", help="CIF file to parse.")
    parser.add_argument(
        "-o", "--output",
        help="Directory to output parsed files. Default: input file directory."
    )
    parser.add_argument(
        "-f", "--filter-by",
        default=SymmetryClass.FAMILY,
        type=SymmetryClass,
        help=(
            "How to separate parsed structures. Supports following arguments:\n"
            "- 'family' (default): parse into 6 crystal families "
            "(triclinic, monoclinic, orthorhombic, tetragonal, hexagonal, cubic)\n"
            "- 'system': parse into 7 crystal systems "
            "(hexagonal family is separated between trigonal and hexagonal systems)\n"
            "- 'pointgroup': parse into 32 crystallographic point groups\n"
            "- 'spacegroup': parse into 230 crystallographic space groups"
        )
    )
    parser.add_argument(
        "-sk", "--special-keys",
        nargs="*",
        help=(
            "CIF Labels to store into structure properties, e.g. can be used "
            "to save structure identifiers attached to it throughout its manipulation "
            "as a python object. Pass 'header' to save structures header (after 'data_')"
        )
    )
    parser.add_argument(
        "--workers", "-w", type=int, metavar="int",
        help=(
            "Number of processes to use in parallel. If not given, will use default of "
            "`tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()` "
            "and execute sequentially."
        )
    )
    args = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at {path.basename(__file__)!r} script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats="cif")
    args_dict["input_file"] = path.abspath(path.realpath(args_dict["input_file"]))
    default_output = path.dirname(args_dict["input_file"])
    args_dict.setdefault("output", default_output)
    check_file_or_dir(args_dict["output"], "dir")
    args_dict["output"] = path.abspath(path.realpath(args_dict["output"]))

    args_dict.setdefault("sequential", False)

    if args_dict.get("special_keys") is not None:
        check_type(args_dict.get("special_keys"), "special_keys", (list,))
        for idx, elt in enumerate(args_dict["special_keys"]):
            check_type(elt, f"special_keys[{idx}]", (str,))

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">", 0)

    return args_dict


def _print_config(args_dict: dict[str, tp.Any]) -> None:
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
    print(f"FILTER CLASS: {args_dict['filter_by'].value!r}")
    print(f"IS SEQUENTIAL: {args_dict.get('sequential')}")
    print(
        "NUMBER OF WORKERS: "
        f"{'auto' if args_dict.get('workers') is None else args_dict.get('workers')} "
        f"{'(ignored)' if args_dict.get('sequential') else ''}"
    )
    print(" ")
    print("------------------------------")
    print(" ")


def get_sym_and_refined_struct(
    structure: Structure, sym_class: SymmetryClass
) -> tuple[str, Structure]:
    """Get the symmetry class and conventional refined version of the structure."""
    spga = SpacegroupAnalyzer(structure)
    refined_struct = spga.get_refined_structure()
    refined_struct.properties = structure.properties

    match sym_class:
        case SymmetryClass.FAMILY:
            family = (
                "hexagonal" if spga.get_crystal_system() == "trigonal"
                else str(spga.get_crystal_system())
            )
            return family, refined_struct

        case SymmetryClass.SYSTEM:
            system = str(spga.get_crystal_system())
            return system, refined_struct

        case SymmetryClass.POINTGROUP:
            pg = spga.get_point_group_symbol().strip()
            return pg, refined_struct

        case SymmetryClass.SPACEGROUP:
            spg = str(spga.get_space_group_number())
            return spg, refined_struct

        case _:
            check_type(sym_class, "sym_class", (SymmetryClass,))
            raise NotImplementedError(f"'sym_class' argument value {sym_class!r} is not supported.")


def main(standalone: bool = True, **kwargs) -> None:
    """
    Find symmetry of CIF structures, refine positions according to symmetry
    and store them in different files depending on their symmetry class.

    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    input_file: str | Path
        CIF file to parse.

    output: str | Path
        Directory to output parsed files. Default: input file directory.

    filter_by: 'family', 'system', 'pointgroup', 'spacegroup'
        How to separate parsed structures. Supports following arguments:
        - 'family' (default): parse into 6 crystal families
        (triclinic, monoclinic, orthorhombic, tetragonal, hexagonal, cubic);
        - 'system': parse into 7 crystal systems
        (hexagonal family is separated between trigonal and hexagonal systems);
        - 'pointgroup': parse into 32 crystallographic point groups;
        - 'spacegroup': parse into 230 crystallographic space groups.

    special_keys: list[str], optional
        CIF Labels to store into structure properties, e.g. can be used to save
        structure identifiers attached to it throughout its manipulation as a python
        object.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.
    """
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)
    _print_config(args)

    cif_file = CIFFile.from_file(
        args["input_file"],
        special_keys=args.get("special_keys"),
        workers=args.get("workers")
    )
    structs, _ = cif_file.parse_structures()

    partial_fn = ft.partial(get_sym_and_refined_struct, sym_class=args["filter_by"])
    sym_desc = "Computing symmetries"
    if args["sequential"]:
        structs = VisualIterator(structs, desc=sym_desc, percent=True, unit="structs")
        structs_tups = [partial_fn(struct) for struct in structs]

    else:
        structs_tups = list(
            process_map(
                partial_fn,
                structs,
                max_workers=args.get("workers"),
                chunksize=min(10, len(structs) // 100 + 1),
                desc=sym_desc
            )
        )

    match args["filter_by"]:
        case SymmetryClass.FAMILY:
            classes = (
                "triclinic", "monoclinic", "orthorhombic", "tetragonal",
                "hexagonal", "cubic"
            )
        case SymmetryClass.SYSTEM:
            classes = (
                "triclinic", "monoclinic", "orthorhombic", "tetragonal",
                "trigonal", "hexagonal", "cubic"
            )
        case SymmetryClass.POINTGROUP:
            classes = tuple(PG_TO_SYSTEM.keys())

        case SymmetryClass.SPACEGROUP:
            classes = tuple(map(str, range(1, 231)))

        case _:
            check_type(args["filter_by"], "--filter-by", (SymmetryClass,))
            raise NotImplementedError(
                f"'--filter-by' argument value {args['filter_by']} is not supported."
            )

    cif_file.clear()
    for sym_class in VisualIterator(
        classes, desc="Parsing symmetry classes", percent=True, unit="symmetry classes"
    ):
        sym_structs = list(filter(lambda t: t[0] == sym_class, structs_tups))
        outfile = path.join(
            args["output"],
            str(path.basename(args["input_file"])).replace(".cif", f"_{sym_class}.cif")
        )
        cif_file.add_structures(
            [struct for _, struct in sym_structs],
            symmetrize=False
        )
        cif_file.write_file(outfile)
        cif_file.clear()
        structs_tups = list(filter(lambda t: t[0] != sym_class, structs_tups))

    assert len(structs_tups) == 0, (
        "Not all structures were parsed at the end, "
        "some symmetry classes may not have been accounted for in the code. "
        f"Unparsed classes: {', '.join(sorted(set(tup[0] for tup in structs_tups)))}."
    )
    print("All structures were properly parsed !")

if __name__ == "__main__":
    main()