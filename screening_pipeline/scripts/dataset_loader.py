#!/usr/bin/python
"""Get MP and/or process OQMD database entries for the pipeline."""

import argparse
import os

# LOCAL IMPORTS
from screening_pipeline.utils import (
    check_file_or_dir, mp_api_download, process_oqmd_json_file
)


########################################


def main() -> None:
    """Main entry point."""

    # ARGUMENTS PARSING BLOCK

    parser = argparse.ArgumentParser(
        description="Get MP and/or process OQMD database entries for the pipeline."
    )

    parser.add_argument(
        "--from-mp-api",
        help=(
            "Path to the JSON file to create and store the Materials Project data in."
            "Fetches entries with thermodynamic data via the Materials Project's REST API, "
            "and store them in a json file usable for the phase_diag_energies.py script at the "
            "--reference argument."
        ),
        metavar="path"
    )
    parser.add_argument(
        "-k", "--mp-api-key",
        help=(
            "If the --from-mp-api arg is used, you can either enter manually a valid MP "
            "API key here or set 'PMG_MAPI_KEY' in .pmgrc.yaml for the program to be able "
            "to fetch the data from the Materials Project API."
        )
    )
    parser.add_argument(
        "--process-oqmd",
        help=(
            "Path to JSON file in OQMD RESTful API format with at least the 'entry_id', "
            "'composition', 'natoms' and 'delta_e' fields. The data will be processed to "
            "fit data formatting needed by phase stability computing script. "
            "The processed data is written to a file with the same name as the one given "
            "but with a '_proc' suffix, at the same location as the one given."
        ),
        metavar="path"
    )

    args = parser.parse_args()


    # MAIN BLOCK

    if args.from_mp_api is not None:
        dir_path = os.path.dirname(args.from_mp_api)
        check_file_or_dir(dir_path, "dir")
        mp_api_download(args.from_mp_api, api_key=args.mp_api_key)

    if args.process_oqmd is not None:
        check_file_or_dir(args.process_oqmd, "file", allowed_formats="json")
        process_oqmd_json_file(path=args.process_oqmd)


if __name__ == "__main__":
    main()
