#!/usr/bin/python
import argparse


def main():
    parser = argparse.ArgumentParser(
        "A command line tool to calculate S.U.N. and RMSD metrics"
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
    parser.add_argument(
        "-t",
        "--threshold",
        default=0.9,
        type=float,
    )

    args = parser.parse_args()

    import os
    import json

    import numpy as np
    from pymatgen.analysis.structure_matcher import StructureMatcher

    from screening_pipeline.utils.vasp_io import batch_extract_vasp_structures
    from screening_pipeline.utils.cif_io import read_cif, extract_cif_from_file
    from screening_pipeline.utils.matcher import remove_equivalent
    from screening_pipeline.utils.ml_vectors import vectors_from_alignn
    from screening_pipeline.utils.distribution import (
        recall,
        precision,
        frechet_distance,
        wasserstein_distance,
    )
    from screening_pipeline.utils.density import get_densities
    from screening_pipeline.utils.rmsd import rmsd_from_structures

    generated, _, _ = read_cif(
        args.generated, keep_rare_gases=True, keep_rare_earths=True
    )
    symmetrized, _, _ = read_cif(
        args.symmetrized, keep_rare_gases=True, keep_rare_earths=True
    )
    dataset, _, _ = read_cif(args.dataset, keep_rare_gases=True, keep_rare_earths=True)

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
    num_novel_unique_stable = sum(map(lambda x: x["stable"], summary))

    # RMSD
    in_structs = [s for s, _ in vasp_structures]
    out_structs = [s for _, s in vasp_structures]
    rmsd = np.mean(rmsd_from_structures(in_structs, out_structs))

    # machine learning
    latent_dataset = vectors_from_alignn(dataset,output="latent")
    latent_gen = vectors_from_alignn(generated,output="latent")

    p = precision(latent_gen, latent_dataset, args.threshold)
    r = recall(latent_gen, latent_dataset, args.threshold)
    fd = frechet_distance(latent_gen, latent_dataset)

    energy_dataset = vectors_from_alignn(dataset,output="energy")
    energy_gen = vectors_from_alignn(generated,output="energy")
    emd_energy = wasserstein_distance(energy_dataset,energy_gen)

    densities_dataset=get_densities(dataset)
    densities_generated=get_densities(generated)
    emd_density = wasserstein_distance(densities_dataset,densities_generated)


    metrics = {
        "dft": {
            "num_generated": num_generated,
            "num_unique": num_unique,
            "num_novel": num_novel,
            "num_novel_unique": num_novel_unique,
            "num_novel_unique_stable": num_novel_unique_stable,
            "SUN": num_novel_unique_stable / num_generated,
            "RMSD": rmsd,
        },
        "ml": {
            "precision": p,
            "recall": r,
            "frechet_distance": fd,
            "EMD_energy":emd_energy,
            "EMD_density":emd_density
        },
    }

    with open(args.output, "w") as fp:
        json.dump(metrics, fp, indent=4)

    print(json.dumps(metrics, indent=4))


if __name__ == "__main__":
    main()
