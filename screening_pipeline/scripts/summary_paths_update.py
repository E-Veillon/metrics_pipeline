#!/usr/bin/python

import argparse

def main() -> None:
    """Main function."""
    
    # ARGUMENTS PARSING BLOCK

    parser = argparse.ArgumentParser(
        description="A command-line tool to modify path designation in JSON summary files."
    )
    parser.add_argument(
        "input_file", help="Path to the JSON file to update paths in."
    )
    parser.add_argument(
        "new_path", help="Path to the directory where structure directories are actually stored."
    )
    parser.add_argument(
        "-a", "--absolute",
        action="store_true",
        help="The new path is abolutized before replacing the old one."
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to spawn."
    )

    args = parser.parse_args()

    import os
    import json
    from typing import List, Dict
    from tqdm.contrib.concurrent import process_map

    assert args.input_file.endswith(".json"), (
        f"'input_file' must be a JSON file."
    )
    assert os.path.isdir(args.new_path), (
        f"'{args.new_path}': No such directory found."
    )
    assert args.workers >= 1, (
        f"'workers' arg must be strictly positive."
    )

    # MAIN BLOCK

    with open(args.input_file, "rt") as fp:
        data = json.load(fp)
    
    assert (
        isinstance(data, List)
        and all([isinstance(struct, Dict) for struct in data])
    ), "The data inside the JSON file must be a list of structure dicts."

    if args.absolute:
        new_path = os.path.abspath(args.new_path)
    else:
        new_path = args.new_path
    
    for struct in data:
        old_path = struct["path"]
        struct_name = os.path.basename(old_path)
        struct["path"] = os.path.join(new_path, struct_name)
    
    with open(args.input_file, "wt") as fp:
        json.dump(data, fp, indent=4)
    
    print(f"the file '{os.path.basename(args.input_file)}' was successfully modified.")

if __name__ == "__main__":
    main()