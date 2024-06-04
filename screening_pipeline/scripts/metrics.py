#!/usr/bin/python
import argparse

from screening_pipeline.utils.crystalnn import to_crystalnn_fingerprint


def main():
    parser = argparse.ArgumentParser(
        "A command line tool to compute S.U.N., RMSD, Coverage Recall and Precision "
        "(COV-R, COV-P), and Earth Mover's Distance (EMD) on densities and energies.\n"
        "S.U.N. metrics require --dataset, --generated, --summary and --preprocessed args.\n"
        "Coverage and EMD metrics require --dataset and --generated args.\n"
        "RMSD metric only requires --summary arg.\n"
        "If some args are not given, corresponding metrics computations will be skipped."
    )

    parser.add_argument(
        "-d", "--dataset", help="Cif file containing the list of known structures."
    )
    parser.add_argument(
        "-g", "--generated", help="Cif file containing all the generated structures."
    )
    parser.add_argument(
        "--no-rare-gas-check",
        action="store_true",
        help=(
            "A flag to disable elimination of structures containing rare gas elements "
            "in the generated structures set.\n"
            "This flag should only be used if it was also used for the preprocessing step."
        ),
        dest="no_rare_gas_check",
    )
    parser.add_argument(
        "--no-rare-earth-check",
        action="store_true",
        help=(
            "A flag to disable elimination of structures containing f-block elements "
            "in the generated structures set.\n"
            "This flag should only be used if it was also used for the preprocessing step."
        ),
        dest="no_rare_earth_check",
    )
    parser.add_argument(
        "-p", "--preprocessed", help="Cif file containing the preprocessed structures."
    )
    parser.add_argument(
        "-s",
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
        help="Threshold for the computation of coverage recall and coverage precision metrics.",
    )

    args = parser.parse_args()

    import os
    import json
    from copy import deepcopy

    import numpy as np
    from pymatgen.io.vasp.inputs import Poscar

    from screening_pipeline.utils import (
        batch_extract_vasp_structures,
        converged_Vasprun,
        read_cif,
        batch_group_by_equivalence,
        remove_equivalent,
        batch_get_novel_structures,
        vectors_from_alignn,
        recall,
        precision,
        frechet_distance,
        wasserstein_distance,
        get_densities,
        rmsd_from_structures,
    )

    if args.dataset is not None:
        print("Loading test set...")
        dataset, _, _ = read_cif(
            filename=args.dataset,
            workers=args.workers,
            keep_rare_gases=True,
            keep_rare_earths=True,
        )
        print("Test set loaded.")
        # remove duplicate structures from the dataset
        dataset, _ = remove_equivalent(
            structures=dataset, workers=args.workers, keep_equivalent=False
        )

    if args.generated is not None:
        print("Loading generated structures...")
        full_generated, _, _ = read_cif(
            filename=args.generated,
            workers=args.workers,
            keep_rare_gases=True,
            keep_rare_earths=True,
        )
        print("Generated structures loaded.")
        if not args.no_rare_gas_check or not args.no_rare_earth_check:
            print("Pruning undesirable elements from generated structures...")
            pruned_generated, _, _ = read_cif(
                filename=args.generated,
                workers=args.workers,
                keep_rare_gases=args.no_rare_gas_check,
                keep_rare_earths=args.no_rare_earth_check,
            )
            print("Pruning finished.")
        else:
            pruned_generated = deepcopy(full_generated)

    if args.preprocessed is not None:
        print("Loading preprocessed structures...")
        preprocessed, _, _ = read_cif(
            filename=args.preprocessed,
            workers=args.workers,
            keep_rare_gases=True,
            keep_rare_earths=True,
        )
        print("Preprocessed structures loaded.")

    if args.summary is not None:
        print("Loading summary file...")
        with open(args.summary, "r") as fp:
            summary = json.load(fp)
        print("Summary file loaded.")

        print("Convert data from summary file to structures...")
        vasp_structures = batch_extract_vasp_structures(
            calc_dirs=[struct["path"] for struct in summary], workers=args.workers
        )
        print("Data converted.")

    dft_metrics = dict.fromkeys(
        (
            "num_generated",
            "num_generated_wo_rare",
            "num_unique", "num_novel",
            "num_unique_novel", "percent_unique_novel",
            "num_stable", "percent_stable",
            "SUN",
            "RMSD"
        )
    )

    ml_metrics = dict.fromkeys(
        ("precision", "recall", "frechet_distance", "EMD_energy", "EMD_density")
    )

    if (
        args.dataset is not None
        and args.generated is not None
        and args.preprocessed is not None
        and args.summary is not None
    ):
        # S.U.N. metrics
        print("Computing S.U.N. metrics...")

        # total count
        dft_metrics["num_generated"] = len(full_generated)
        dft_metrics["num_generated_wo_rare"] = len(pruned_generated)

        # Unique count
        dft_metrics["num_unique"] = len(preprocessed)

        # novel count
        concat_novel, _ = remove_equivalent(
            structures=pruned_generated + dataset, workers=args.workers, keep_equivalent=False
        )
        dft_metrics["num_novel"] = len(concat_novel) - len(dataset)

        # novel + unique count
        concat_novel_unique, _ = remove_equivalent(
            structures=preprocessed + dataset, workers=args.workers, keep_equivalent=False
        )
        dft_metrics["num_unique_novel"] = len(concat_novel_unique) - len(dataset)
        dft_metrics["percent_unique_novel"] = dft_metrics["num_unique_novel"] / len(pruned_generated)

        # stable count
        dft_metrics["num_stable"] = sum(map(lambda x: x["stable"], summary))
        dft_metrics["percent_stable"] = dft_metrics["num_stable"] / len(vasp_structures)

        # S.U.N. percentage
        dft_metrics["SUN"] = dft_metrics["percent_stable"] * dft_metrics["percent_unique_novel"]

        print("S.U.N. metrics computed.")
        for key, val in dft_metrics.items():
            if key == "RMSD":
                continue
            print(f"{key} = {val}")

    if args.summary is not None:
        # RMSD metric
        print("Computing RMSD metric...")
        in_structs = [s for s, _ in vasp_structures]
        out_structs = [s for _, s in vasp_structures]
        dft_metrics["RMSD"] = np.mean(
            rmsd_from_structures(in_structs, out_structs)
        ).item()
        print("RMSD metric computed.")
        print(f"RMSD = {dft_metrics['RMSD']}")

    if args.dataset is not None and args.generated is not None:
        # machine learning metrics (COV-R, COV-P, energy EMD, density EMD)
        print("Computing latent space metrics (COV-R, COV-P)...")
        fingerprint_dataset = to_crystalnn_fingerprint(dataset, workers=args.workers)
        fingerprint_gen = to_crystalnn_fingerprint(full_generated, workers=args.workers)

        fingerprint_dataset, fingerprint_gen = map(np.array,zip(
            *filter(
                lambda x: x[0] is not None and x[1] is not None,
                zip(fingerprint_dataset, fingerprint_gen),
            )
        ))

        ml_metrics["precision"] = precision(
            fingerprint_gen, fingerprint_dataset, args.threshold
        )
        ml_metrics["recall"] = recall(
            fingerprint_gen, fingerprint_dataset, args.threshold
        )

        latent_dataset = vectors_from_alignn(dataset, output="latent")
        latent_gen = vectors_from_alignn(full_generated, output="latent")
        ml_metrics["frechet_distance"] = frechet_distance(latent_gen, latent_dataset)
        print("Latent space metrics computed.")
        print(f"COV-P = {ml_metrics['precision']}")
        print(f"COV-R = {ml_metrics['recall']}")
        print(f"Frechet Distance = {ml_metrics['frechet_distance']}")

        print("Computing properties EMD metrics...")
        energy_dataset = vectors_from_alignn(dataset, output="energy")
        energy_gen = vectors_from_alignn(full_generated, output="energy")
        ml_metrics["EMD_energy"] = wasserstein_distance(energy_dataset, energy_gen)

        densities_dataset = get_densities(dataset)
        densities_generated = get_densities(full_generated)
        ml_metrics["EMD_density"] = wasserstein_distance(
            densities_dataset, densities_generated
        )
        print("Properties EMD metrics computed.")
        print(f"Density EMD = {ml_metrics['EMD_density']}")
        print(f"Energy EMD = {ml_metrics['EMD_energy']}")

    metrics = {"dft": dft_metrics, "ml": ml_metrics}

    print("Writing output file...")

    with open(args.output, "w") as fp:
        json.dump(metrics, fp, indent=4)

    print(f"Output file successfully written at location '{args.output}'.")
    print(json.dumps(metrics, indent=4))


if __name__ == "__main__":
    main()
