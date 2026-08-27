"""Concatenate JSON structure dataset files seemlessly."""

import argparse as ap

from core.genmat_io import StructureFile


def main() -> None:
    parser = ap.ArgumentParser()
    parser.add_argument("files", nargs="+", help="JSON databases to concatenate")
    parser.add_argument(
        "-o", "--output", default="concatenated.json",
        help="JSON output file to write concatenated dataset. Defaults to %(default)s."
    )
    args = parser.parse_args()

    cat_file = StructureFile.from_file(args.files[0])

    for file in args.files[1:]:
        sfile = StructureFile.from_file(file)
        cat_file = cat_file + sfile

    cat_file.write_file(args.output)


if __name__ == "__main__":
    main()

