#!/usr/bin/python
import argparse


def main():
    parser = argparse.ArgumentParser(
        "A command line tool to calculate S.U.N. and RSMD metrics"
    )

    parser.add_argument(
        "-g", "--generated", help="Cif file containing all the generated structures."
    )
    parser.add_argument(
        "-s", "--symmetrized", help="Cif file containing the symetrized structures."
    )
    parser.add_argument(
        "-d", "--dataset", help="Cif file containing the list of known structures."
    )
    parser.add_argument(
        "-S",
        "--summary",
        help="Json file containing the summary generated after the stability screening.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="metrics.json",
        help="Output file containing the calculated metrics (json format).",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to spawn for parallelized steps.",
        metavar="int",
    )

    args = parser.parse_args()

    import os
    import json

    import numpy as np
    from pymatgen.analysis.structure_matcher import StructureMatcher

    from screening_pipeline.utils.vasp_io import batch_extract_vasp_structures
    from screening_pipeline.utils.cif_io import read_cif, extract_cif_from_file
    from screening_pipeline.utils.matcher import remove_equivalent

    generated = read_cif(args.generated, keep_rare_gases=True, keep_rare_earths=True)
    symmetrized = read_cif(
        args.symmetrized, keep_rare_gases=True, keep_rare_earths=True
    )
    dataset = read_cif(args.dataset, keep_rare_gases=True, keep_rare_earths=True)

    with open(args.summary, "r") as fp:
        summary = json.load(fp)

    vasp_structures = batch_extract_vasp_structures(
        [struct["path"] for struct in summary]
    )

    # remove duplicate structures from the dataset
    dataset, _ = remove_equivalent(
        structures=dataset, workers=args.workers, keep_equivalent=False
    )

    # total
    num_generated = len(generated)

    # unique count
    num_unique = len(symmetrized)

    # novel count
    concat_novel, _ = remove_equivalent(
        structures=generated + dataset, workers=args.workers, keep_equivalent=False
    )
    num_novel = len(concat_novel) - len(dataset)

    # novel + unique count
    concat_novel_unique, _ = remove_equivalent(
        structures=symmetrized + dataset, workers=args.workers, keep_equivalent=False
    )
    num_novel_unique = len(concat_novel_unique) - len(dataset)

    # novel + unique + stable count
    num_novel_unique_stable = len(filter(map(lambda x: x["stable"], summary)))

    # RMSD
    rmsd_list = []
    matcher = StructureMatcher()
    for in_struct, out_struct in vasp_structures:
        if in_struct is None or out_struct is None:
            continue

        rms = matcher.get_rms_dist(in_struct, out_struct)

        if rms is not None:
            rmsd_list.append(rms[0])
    rmsd = float(np.mean(rmsd_list))

    metrics = {
        "num_generated": num_generated,
        "num_unique": num_unique,
        "num_novel": num_novel,
        "num_novel_unique": num_novel_unique,
        "num_novel_unique_stable": num_novel_unique_stable,
        "SUN": num_novel_unique_stable / num_generated,
        "rmsd": rmsd,
    }

    with open(args.output, "w") as fp:
        json.dump(metrics, fp, indent=4)

    print(json.dum(metrics, indent=4))


if __name__ == "__main__":
    main()
