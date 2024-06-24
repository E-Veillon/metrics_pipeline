#!/usr/bin/python

import argparse
import json
import numpy as np

# LOCAL IMPORTS
from screening_pipeline.utils import (
    check_file_format, check_file_or_dir,
    batch_extract_vasp_structures, read_cif,
    remove_equivalent, vectors_from_alignn,
    recall, precision, frechet_distance,
    emd_wrapper, get_densities,
    rmsd_from_structures, to_crystalnn_fingerprint
)


########################################


def _assert_args(args: argparse.Namespace) -> None:
    """Check arguments values validity."""
    if args.dataset is not None:
        check_file_or_dir(args.dataset, "file", format="cif")
    if args.generated is not None:
        check_file_or_dir(args.generated, "file", format="cif")
    if args.uniques is not None:
        check_file_or_dir(args.uniques, "file", format="cif")
    if args.valid is not None:
        check_file_or_dir(args.valid, "file", format="cif")
    if args.summary is not None:
        check_file_or_dir(args.summary, "file", format="json")

    check_file_format(args.output, format="json")

    if args.workers < 1:
        raise ValueError(
            f"'workers' argument must be strictly positive (got {args.workers})."
        )
    if args.threshold <= 0.0:
        raise ValueError(
            f"'threshold' argument must be strictly positive (got {args.threshold})."
        )


########################################


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description=
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
        "-u", "--uniques", help="Cif file containing the preprocessed unique structures."
    )
    parser.add_argument(
        "-v", "--valid", help="Cif file containing the preprocessed valid structures."
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
        "--test-min-vol",
        action="store_true",
        help=(
            "A debug flag to assume unicity of unlikely structures having a volume under "
            "1 Angström^3 without passing them into structure matching, which could cause "
            "the program to be softlocked during S.U.N. computations. Only pass it if such "
            "problems were to arise."
        ),
    )
    parser.add_argument(
        "-t",
        "--threshold",
        default=0.4,
        type=float,
        help="Threshold for the computation of coverage recall and coverage precision metrics.",
    )

    args = parser.parse_args()

    _assert_args(args)

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
        dataset, _, _ = remove_equivalent(
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
        #if not args.no_rare_gas_check or not args.no_rare_earth_check:
        #    print("Pruning undesirable elements from generated structures...")
        #    pruned_generated, _, _ = read_cif(
        #        filename=args.generated,
        #        workers=args.workers,
        #        keep_rare_gases=args.no_rare_gas_check,
        #        keep_rare_earths=args.no_rare_earth_check,
        #    )
        #    print("Pruning finished.")
        #else:
        #    pruned_generated = deepcopy(full_generated)

    if args.uniques is not None:
        print("Loading preprocessed uniques structures...")
        uniques, _, _ = read_cif(
            filename=args.uniques,
            workers=args.workers,
            keep_rare_gases=True,
            keep_rare_earths=True,
        )
        print("Preprocessed uniques structures loaded.")

    if args.valid is not None:
        print("Loading preprocessed valid structures...")
        valids, _, _ = read_cif(
            filename=args.valid,
            workers=args.workers,
            keep_rare_gases=True,
            keep_rare_earths=True,
        )
        print("Preprocessed valid structures loaded.")

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
            "num_generated",#"num_generated_wo_rare",
            "num_valid", "percent_valid",
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

    if args.valid is not None:
        dft_metrics["num_valid"] = len(valids)
        prop_valid = len(valids) / len(full_generated)
        dft_metrics["percent_valid"] = round(prop_valid * 100, 6)

    if (
        args.dataset is not None
        and args.generated is not None
        and args.uniques is not None
        and args.summary is not None
    ):
        # S.U.N. metrics
        print("Computing S.U.N. metrics...")

        # total count
        dft_metrics["num_generated"] = len(full_generated)
        #dft_metrics["num_generated_wo_rare"] = len(pruned_generated)

        # Unique count
        dft_metrics["num_unique"] = len(uniques)

        # novel count
        concat_novel, _, nbr_unmatched = remove_equivalent(
            structures=full_generated + dataset,
            workers=args.workers,
            test_volume=args.test_min_vol,
            keep_equivalent=False
        )
        dft_metrics["num_novel"] = len(concat_novel) - len(dataset)

        if nbr_unmatched != 0:
            dft_metrics["unmatched_novel"] = nbr_unmatched

        # novel + unique count
        concat_novel_unique, _, nbr_unmatched = remove_equivalent(
            structures=uniques + dataset,
            workers=args.workers,
            test_volume=args.test_min_vol,
            keep_equivalent=False
        )
        dft_metrics["num_unique_novel"] = len(concat_novel_unique) - len(dataset)
        prop_unique_novel = dft_metrics["num_unique_novel"] / len(full_generated)
        dft_metrics["percent_unique_novel"] = round(prop_unique_novel * 100, 6)

        if nbr_unmatched != 0:
            dft_metrics["unmatched_novel_unique"] = nbr_unmatched

        # stable count
        dft_metrics["num_stable"] = sum(map(lambda x: x["stable"], summary))
        prop_stable = dft_metrics["num_stable"] / len(full_generated)
        dft_metrics["percent_stable"] = round(prop_stable * 100, 6)

        # S.U.N. percentage
        prop_sun = prop_stable * prop_unique_novel
        dft_metrics["S.U.N."] = round(prop_sun * 100, 6)

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
        ml_metrics["EMD_energy"] = emd_wrapper(energy_dataset, energy_gen)

        densities_dataset = get_densities(dataset)
        densities_generated = get_densities(full_generated)
        ml_metrics["EMD_density"] = emd_wrapper(
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
