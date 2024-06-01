#!/usr/bin/python
import argparse


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
        dest="no_rare_gas_check"
    )
    parser.add_argument(
        "--no-rare-earth-check",
        action="store_true",
        help=(
            "A flag to disable elimination of structures containing f-block elements "
            "in the generated structures set.\n"
            "This flag should only be used if it was also used for the preprocessing step."
        ),
        dest="no_rare_earth_check"
    )
    parser.add_argument(
        "-p", "--preprocessed", help="Cif file containing the preprocessed structures."
    )
    parser.add_argument(
        "-s", "--summary",
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
        help="Threshold for the computation of coverage recall and coverage precision metrics."
    )

    args = parser.parse_args()

    import os
    import json
    from copy import deepcopy

    import numpy as np
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.analysis.structure_matcher import StructureMatcher

    from screening_pipeline.utils.vasp_io import batch_extract_vasp_structures, converged_Vasprun
    from screening_pipeline.utils.cif_io import read_cif
    from screening_pipeline.utils.matcher import batch_group_by_equivalence, remove_equivalent, flatten
    from screening_pipeline.utils.ml_vectors import vectors_from_alignn
    from screening_pipeline.utils.distribution import (
        recall,
        precision,
        frechet_distance,
        wasserstein_distance,
    )
    from screening_pipeline.utils.density import get_densities
    from screening_pipeline.utils.rmsd import rmsd_from_structures

    if args.dataset is not None:
        print("Loading test set...")
        dataset, _, _ = read_cif(
            filename=args.dataset,
            workers=args.workers,
            keep_rare_gases=True,
            keep_rare_earths=True
        )
        print("Test set loaded.")
        # remove duplicate structures from the dataset
        dataset, _ = remove_equivalent(
        structures=dataset, workers=args.workers, keep_equivalent=False
        )

    if args.generated is not None:
        print("Loading genrated structures...")
        full_generated, _, _ = read_cif(
            filename=args.generated,
            workers=args.workers,
            keep_rare_gases=True,
            keep_rare_earths=True
        )
        print("Generated structures loaded.")
        if not args.no_rare_gas_check or not args.no_rare_earth_check:
            print("Pruning undesirable elements from generated structures...")
            pruned_generated, _, _ = read_cif(
                filename=args.generated,
                workers=args.workers,
                keep_rare_gases=args.no_rare_gas_check,
                keep_rare_earths=args.no_rare_earth_check
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
            keep_rare_earths=True
        )
        print("Preprocessed structures loaded.")

    if args.summary is not None:
        print("Loading summary file...")
        with open(args.summary, "r") as fp:
            summary = json.load(fp)
        print("Summary file loaded.")

        print("Convert data from summary file to structures...")
        vasp_structures = batch_extract_vasp_structures(
            calc_dirs=[struct["path"] for struct in summary],
            workers=args.workers
        )
        print("Data converted.")

    dft_metrics = dict.fromkeys(
        (
            "num_generated", "num_generated_wo_rare", "num_unique",
            #"num_novel", "num_novel_unique",
            "num_stable_unique", "num_stable_unique_novel",
            "SUN", "RMSD"
        )
    )

    ml_metrics = dict.fromkeys(
        (
            "precision", "recall", "frechet_distance",
            "EMD_energy", "EMD_density"
        )
    )

    if (args.dataset is not None and args.generated is not None and
        args.preprocessed is not None and args.summary is not None):
        # S.U.N. metrics
        print("Computing S.U.N. metrics...")

        # total count
        dft_metrics["num_generated"] = len(full_generated)
        dft_metrics["num_generated_wo_rare"] = len(pruned_generated)

        # Unique count
        dft_metrics["num_unique"] = len(preprocessed)

        """
        # novel count
        concat_novel, _ = remove_equivalent(
            structures=pruned_generated + dataset, workers=args.workers, keep_equivalent=False
        )
        dft_metrics["num_novel"] = len(concat_novel) - len(dataset)

        # novel + unique count
        concat_novel_unique, _ = remove_equivalent(
            structures=symmetrized + dataset, workers=args.workers, keep_equivalent=False
        )
        dft_metrics["num_novel_unique"] = len(concat_novel_unique) - len(dataset)

        # novel + unique + stable count
        num_novel_unique_stable = sum(map(lambda x: x["stable"], summary))
        """

        # Stable + Unique count
        paths_stable_unique = list(filter(lambda data: data["stable"].lower() == "true", summary))
        dft_metrics["num_stable_unique"] = len(paths_stable_unique)

        # Stable + Unique + Novel count (S.U.N.)
        structs_stable_unique = list(
            map(
                lambda data: converged_Vasprun(data["path"]).initial_structure,
                paths_stable_unique
            )
        )

        concat_sun = batch_group_by_equivalence(
            structures=structs_stable_unique + dataset,
            workers=args.workers,
            comment="Comparing known and generated structures"
        )

        sun_structs = []

        for comp_group in concat_sun:
            sun_structs += list(filter(
                lambda l: len(l) == 1 and l[0] not in dataset,
                comp_group
            ))
        sun_structs = flatten(sun_structs)

        dft_metrics["num_stable_unique_novel"] = len(sun_structs)
        print("S.U.N. metrics computed.")

    if args.summary is not None:
        # RMSD metric
        print("Computing RMSD metric...")
        in_structs = [s for s, _ in vasp_structures]
        out_structs = [s for _, s in vasp_structures]
        dft_metrics["rmsd"] = np.mean(rmsd_from_structures(in_structs, out_structs))
        print("RMSD metric computed.")

    if args.dataset is not None and args.generated is not None:
        # machine learning metrics (COV-R, COV-P, energy EMD, density EMD)
        print("Computing latent space metrics (COV-R, COV-P)...")
        latent_dataset = vectors_from_alignn(dataset,output="latent")
        latent_gen = vectors_from_alignn(full_generated,output="latent")

        ml_metrics["precision"] = precision(latent_gen, latent_dataset, args.threshold)
        ml_metrics["recall"] = recall(latent_gen, latent_dataset, args.threshold)
        ml_metrics["frechet_distance"] = frechet_distance(latent_gen, latent_dataset)
        print("Latent space metrics computed.")

        print("Computing properties EMD metrics...")
        energy_dataset = vectors_from_alignn(dataset,output="energy")
        energy_gen = vectors_from_alignn(full_generated,output="energy")
        ml_metrics["EMD_energy"] = wasserstein_distance(energy_dataset,energy_gen)

        densities_dataset=get_densities(dataset)
        densities_generated=get_densities(full_generated)
        ml_metrics["EMD_density"] = wasserstein_distance(densities_dataset,densities_generated)
        print("Properties EMD metrics computed.")

    metrics = {"dft": dft_metrics, "ml": ml_metrics}

    print("Writing output file...")

    with open(args.output, "w") as fp:
        json.dump(metrics, fp, indent=4)

    print(f"Output file successfully written at location '{args.output}'.")
    print(json.dumps(metrics, indent=4))


if __name__ == "__main__":
    main()
