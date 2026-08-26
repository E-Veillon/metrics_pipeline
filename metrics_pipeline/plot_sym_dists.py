"""Generate histograms from JSON symmetry mapping files produced by the symmetry_mapping.py script."""

import argparse as ap
import matplotlib.pyplot as plt

from src.genmat_io import JsonLoader, check_file_or_dir, check_file_format


def main() -> None:
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", help="JSON file to plot.")
    parser.add_argument("-o", "--output", help="Path to save the plot as PNG file.")
    parser.add_argument(
        "-c", "--category", default="system", choices=["system", "point_group", "space_group"],
        help="Level of symmetry to consider in the plot."
    )
    parser.add_argument(
        "-t", "--total", action=ap.BooleanOptionalAction, default=False,
        help="Add a 'total' bar containing the sum of all categories."
    )
    args = parser.parse_args()

    if args.output is None:
        args.output = args.input_file.replace(".json", f"_{args.category}.png")

    check_file_or_dir(args.input_file, "file", allowed_formats="json")
    check_file_format(args.output, allowed_formats="png")

    data = JsonLoader(args.input_file).load_as_dict()

    level = f"by_{args.category}"
    name = args.category.replace("_", " ")
    subdata = [(system, total["total"]) for system, total in data[level].items()]
    categories = [t[0] for t in subdata] + ["uncomputable"]
    values = [t[1] for t in subdata] + [data["uncomputable"]["total"]]

    if args.total:
        categories.append("total")
        values.append(sum(values))

    plt.figure()

    bars = plt.bar(categories, values)
    plt.bar_label(bars)

    plt.xlabel(name)
    plt.ylabel("number")
    plt.title(f"Crystal {name} ({args.input_file})")
    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig(args.output, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
