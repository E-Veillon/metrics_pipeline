"""Concatenate JSON structure dataset files seemlessly."""

import argparse as ap

from src.io import StructureFile


def main() -> None:
    parser = ap.ArgumentParser()
    parser.add_argument("files", nargs="+", help="JSON databases to concatenate")
    parser.add_argument(
        "-o", "--output", default="concatenated.json",
        help="JSON output file to write concatenated dataset. Defaults to %(default)s."
    )
    args = parser.parse_args()

    cat_file = None

    for file in args.files:
        sfile = StructureFile.from_file(file)
        if cat_file is None:
            cat_file = sfile
        else:
            cat_file = cat_file + sfile

    cat_file.write_file(args.output)


if __name__ == "__main__":
    main()

