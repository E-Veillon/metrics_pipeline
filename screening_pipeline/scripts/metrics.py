#!/usr/bin/python
"""
A command line tool to compute various metrics based on previously done computations.
Computable metrics are:\n
- Validity,
- Stability, Unicity, Novelty (S.U.N),\n
- Average Root Mean Square Displacement (RMSD),\n
- Coverage - Precision (COV-P) and Coverage - Recall (COV-R),\n
- Fréchet ALIGNN Distance (FAD),\n
- Earth Mover's Distance (EMD) on energy and density distributions.\n
"""

import os
import json
import argparse
from typing import Optional
import numpy as np

# LOCAL IMPORTS
from screening_pipeline.utils import (
    check_num_value, PathLike, CONFIGPATH,
    check_file_format, check_file_or_dir, yaml_loader,
    batch_extract_vasp_structures, read_cif,
    batch_group_by_equivalence, batch_get_novel_structures,
    remove_equivalent, vectors_from_alignn,
    recall, precision, frechet_distance,
    get_emd, get_densities,
    rmsd_from_structures, to_crystalnn_fingerprint
)


########################################


def _assert_args(args: argparse.Namespace) -> None:
    """Check arguments values validity."""
    config_file = os.path.join(CONFIGPATH, args.config)
    check_file_or_dir(config_file, "file", allowed_formats="yaml")
    if args.dataset is not None:
        check_file_or_dir(args.dataset, "file", allowed_formats="cif")
    if args.generated is not None:
        check_file_or_dir(args.generated, "file", allowed_formats="cif")
    if args.uniques is not None:
        check_file_or_dir(args.uniques, "file", allowed_formats="cif")
    if args.valid is not None:
        check_file_or_dir(args.valid, "file", allowed_formats="cif")
    if args.sun_summary is not None:
        check_file_or_dir(args.sun_summary, "file", allowed_formats="json")
    if args.relax_summary is not None:
        check_file_or_dir(args.relax_summary, "file", allowed_formats="json")
    check_file_format(args.output, allowed_formats="json")
    check_num_value(args.workers, "--workers", ">", 0)
    check_num_value(args.threshold, "--threshold", ">", 0.0)

def _match_file_arg_need(
    arg_name: str, filename: Optional[PathLike] = None, is_needed: bool = False
) -> bool:
    match (filename is not None, is_needed):
        case (False, False):
            return False

        case (True, False):
            print(
                f"--{arg_name} arg was given but no activated metric needs it, "
                "therefore its loading is skipped for efficiency."
            )
            return False

        case (False, True):
            raise ValueError(
                f"Config file activates some metrics that need a {arg_name} file, "
                f"but --{arg_name} arg was not given."
            )

        case (True, True):
            return True

def _print_metrics_config(config: dict) -> None:
    print("- ACTIVATED METRICS -")
    print(" ")
    print(f"Validity: {config.get('Validity')}")
    print(f"Stability, Unicity, Novelty (SUN): {config.get('SUN')}")
    print(f"Avg. Root Mean Square Displacement (RMSD): {config.get('RMSD')}")
    print(f"Coverage - Precision (COV-P): {config.get('COV-P')}")
    print(f"Coverage - Recall (COV-R): {config.get('COV-R')}")
    print(f"Fréchet ALIGNN Distance (FAD): {config.get('FAD')}")
    print(f"Earth Mover's Distance (EMD) on energy: {config.get('EMD_energy')}")
    print(f"Earth Mover's Distance (EMD) on density: {config.get('EMD_density')}")
    print(" ")
    print("------------------------------")
    print(" ")

########################################


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description=
        "A command line tool to compute various metrics based on previously done "
        "computations. Computable metrics are:\n"
        "- Validity,"
        "- Stability, Unicity, Novelty (S.U.N),\n"
        "- Average Root Mean Square Displacement (RMSD),\n"
        "- Coverage - Precision (COV-P) and Coverage - Recall (COV-R),\n"
        "- Fréchet ALIGNN Distance (FAD),\n"
        "- Earth Mover's Distance (EMD) on energy and density distributions.\n"
    )
    parser.add_argument(
        "generated",
        help=(
            "Cif file containing all the generated structures. "
        )
    )
    parser.add_argument(
        "-c", "--config",
        default="metrics_defaults.yaml",
        help=(
            "Configuration file in Yaml format stating which metrics should be computed. "
            f"The file needs to be in {CONFIGPATH} to be found. "
            "Defaults to %(default)s. This argument adds flexibility by allowing to copy "
            "the default config file to make several computations presets as needed and "
            "provide the one needed here."
        )
    )
    parser.add_argument(
        "-d", "--dataset", help=(
            "Cif file containing the list of known structures. "
            "Necessary for: S.U.N., COV-P, COV-R, FAD, EMD(density), EMD(energy)."
        )
    )
    parser.add_argument(
        "-v", "--valid",
        help=(
            "Cif file containing the preprocessed valid structures. "
            "Necessary for: Validity."
        )
    )
    parser.add_argument(
        "-u", "--uniques",
        help=(
            "Cif file containing the preprocessed unique structures."
            "Necessary for: S.U.N."
        )
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
        "-s", "--sun-summary",
        help=(
            "Json file containing the summary of phase diagrams instability energies. "
            "Necessary for: S.U.N."
        )
    )
    parser.add_argument(
        "-r", "--relax-summary",
        help=(
            "Json file containing a summary for the relaxation step. "
            "Necessary for: RMSD."
        )
    )
    parser.add_argument(
        "-o", "--output",
        default="metrics.json",
        help="Output file containing the calculated metrics (json format).",
    )
    parser.add_argument(
        "-w", "--workers",
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
        "-t", "--threshold",
        default=0.4,
        type=float,
        help="Threshold for the computation of coverage recall and coverage precision metrics.",
    )

    args = parser.parse_args()

    _assert_args(args)

    print("===== LOAD NECESSARY DATA FILES =====")

    CONFIG = yaml_loader(os.path.join(CONFIGPATH, args.config), on_error='raise')
    dataset_needed = (
        CONFIG.get("SUN")
        or CONFIG.get("COV-P")
        or CONFIG.get("COV-R")
        or CONFIG.get("FAD")
        or CONFIG.get("EMD_energy")
        or CONFIG.get("EMD_density")
    )
    valid_needed = CONFIG.get("Validity")
    uniques_needed = sun_summary_needed = CONFIG.get("SUN")
    relax_summary_needed = CONFIG.get("RMSD")

    _print_metrics_config(CONFIG)

    print("Loading generated structures...")
    generated, *_ = read_cif(
        filename=args.generated,
        workers=args.workers,
        keep_rare_gases=args.no_rare_gas_check,
        keep_rare_earths=args.no_rare_earth_check
    )
    print("Generated structures loaded.")

    if _match_file_arg_need("dataset", args.dataset, dataset_needed):
        print("Loading dataset...")
        dataset, *_ = read_cif(
            filename=args.dataset,
            workers=args.workers,
            keep_rare_gases=args.no_rare_gas_check,
            keep_rare_earths=args.no_rare_earth_check
        )
        print("Dataset loaded.")
        # remove duplicate structures from the dataset
        dataset, *_ = remove_equivalent(
            structures=dataset, workers=args.workers, keep_equivalent=False
        )

    if _match_file_arg_need("valid", args.valid, valid_needed):
        print("Loading preprocessed valid structures...")
        valids, *_ = read_cif(
            filename=args.valid,
            workers=args.workers,
            keep_rare_gases=args.no_rare_gas_check,
            keep_rare_earths=args.no_rare_earth_check
        )
        print("Preprocessed valid structures loaded.")

    if _match_file_arg_need("uniques", args.uniques, uniques_needed):
        print("Loading preprocessed uniques structures...")
        uniques, *_ = read_cif(
            filename=args.uniques,
            workers=args.workers,
            keep_rare_gases=args.no_rare_gas_check,
            keep_rare_earths=args.no_rare_earth_check
        )
        print("Preprocessed uniques structures loaded.")

    if _match_file_arg_need("sun-summary", args.sun_summary, sun_summary_needed):
        print("Loading stability summary file...")
        with open(args.sun_summary, "rt", encoding="utf-8") as fp:
            sun_summary = json.load(fp)
        print("Stability summary file loaded.")

        print("Convert data from stability summary file to structures...")
        stable_structs_paths = [
            data["path"] for data in filter(lambda d: d["stable"], sun_summary)
        ]
        stable_structs = batch_extract_vasp_structures(
            calc_dirs=stable_structs_paths, workers=args.workers
        )
        stable_structs = [s for s, _ in stable_structs]
        print("Data converted.")

    if _match_file_arg_need("relax-summary", args.relax_summary, relax_summary_needed):
        print("Loading relaxations summary file...")
        with open(args.relax_summary, "rt", encoding="utf-8") as fp:
            relax_summary = json.load(fp)
        print("Relaxations summary file loaded.")

        print("Convert data from relaxations summary file to structures...")
        converged_relax_paths = [
            data["path"] for data in filter(lambda d: d["converged"], relax_summary)
        ]
        relax_structures = batch_extract_vasp_structures(
            calc_dirs=converged_relax_paths,
            workers=args.workers
        )
        print("Data converted.")

    general_metrics = dict.fromkeys(
        ("num_generated", "num_valid", "percent_valid")
    )

    dft_metrics = dict.fromkeys(
        (
            "num_unique", "percent_unique",
            "num_novel", "percent_novel",
            "num_unique_novel", "percent_unique_novel",
            "num_stable", "percent_stable",
            "num_SUN", "percent_SUN",
            "RMSD"
        )
    )

    ml_metrics = dict.fromkeys(
        ("precision", "recall", "frechet_distance", "EMD_energy", "EMD_density")
    )

    print("===== COMPUTE ACTIVATED METRICS =====")

    # total number of generated structures
    general_metrics["num_generated"] = len(generated)

    if CONFIG.get("Validity"):
        # Validity metric
        general_metrics["num_valid"] = len(valids)
        prop_valid = len(valids) / len(generated)
        general_metrics["percent_valid"] = round(prop_valid * 100, 6)

    if CONFIG.get("SUN"):
        # S.U.N. metrics
        print("Computing S.U.N. metrics...")

        # Unique count
        dft_metrics["num_unique"] = len(uniques)
        prop_unique = dft_metrics["num_unique"] / general_metrics["num_generated"]
        dft_metrics["percent_unique"] = round(prop_unique * 100, 6)

        # novel count
        grouped_structs, nbr_unmatched = batch_group_by_equivalence(
            structures=generated + dataset,
            workers=args.workers,
            comment="Compare generated and dataset"
        )
        novel_structs = batch_get_novel_structures(
            grouped_structs, dataset, workers=args.workers
        )
        dft_metrics["num_novel"] = len(novel_structs)
        prop_novel = dft_metrics["num_novel"] / general_metrics["num_generated"]
        dft_metrics["percent_novel"] = round(prop_novel * 100, 6)

        if nbr_unmatched != 0:
            dft_metrics["num_unmatched_novel"] = nbr_unmatched

        # novel + unique count
        grouped_structs, nbr_unmatched = batch_group_by_equivalence(
            structures=uniques + dataset,
            workers=args.workers,
            comment="Compare uniques and dataset"
        )
        unique_novel_structs = batch_get_novel_structures(
            grouped_structs, dataset, workers=args.workers
        )
        dft_metrics["num_unique_novel"] = len(unique_novel_structs)
        prop_unique_novel = dft_metrics["num_unique_novel"] / general_metrics["num_generated"]
        dft_metrics["percent_unique_novel"] = round(prop_unique_novel * 100, 6)

        if nbr_unmatched != 0:
            dft_metrics["num_unmatched_unique_novel"] = nbr_unmatched

        # stable count
        dft_metrics["num_stable"] = len(stable_structs)
        prop_stable = dft_metrics["num_stable"] / general_metrics["num_generated"]
        dft_metrics["percent_stable"] = round(prop_stable * 100, 6)

        # S.U.N. count
        sun_structs = list(
            filter(
                lambda struct:any(
                    struct == uniq_novel for uniq_novel in unique_novel_structs
                ), stable_structs
            )
        )
        dft_metrics["num_SUN"] = len(sun_structs)
        prop_sun = dft_metrics["num_SUN"] / general_metrics["num_generated"]
        dft_metrics["percent_SUN"] = round(prop_sun * 100, 6)

        for key, val in dft_metrics.items():
            if key == "RMSD":
                continue
            print(f"{key} = {val}")

    if CONFIG.get("RMSD"):
        # RMSD metric
        print("Computing RMSD metric...")
        in_structs = [s for s, _ in relax_structures]
        out_structs = [s for _, s in relax_structures]
        dft_metrics["RMSD"] = round(
            np.mean(
                rmsd_from_structures(in_structs, out_structs)
            ).item(),
            ndigits=6
        )
        print(f"RMSD = {dft_metrics['RMSD']}")

    if CONFIG.get("COV-P") or CONFIG.get("COV-R"):
        # Compute Coverage (Precision, Recall)
        print("Computing fingerprints for Coverage metrics...")
        fingerprint_dataset = to_crystalnn_fingerprint(dataset, workers=args.workers)
        fingerprint_gen = to_crystalnn_fingerprint(generated, workers=args.workers)

        fingerprint_dataset, fingerprint_gen = map(np.array,zip(
            *filter(
                lambda x: x[0] is not None and x[1] is not None,
                zip(fingerprint_dataset, fingerprint_gen),
            )
        ))
        if CONFIG.get("COV-P"):
            print("Computing Coverage (Precision)...")
            ml_metrics["precision"] = precision(
                fingerprint_gen, fingerprint_dataset, args.threshold
            )
            print(f"COV-P = {ml_metrics['precision']}")

        if CONFIG.get("COV-R"):
            print("Computing Coverage (Recall)...")
            ml_metrics["recall"] = recall(
                fingerprint_gen, fingerprint_dataset, args.threshold
            )
            print(f"COV-R = {ml_metrics['recall']}")

    if CONFIG.get("FAD"):
        # Compute Fréchet ALIGNN Distance
        print("Computing Fréchet ALIGNN Distance metric...")
        latent_dataset = vectors_from_alignn(dataset, output="latent")
        latent_gen = vectors_from_alignn(generated, output="latent")
        ml_metrics["frechet_distance"] = frechet_distance(latent_gen, latent_dataset)
        print(f"Frechet Distance = {ml_metrics['frechet_distance']}")

    if CONFIG.get("EMD_energy"):
        # Compute Earth Mover's Distance on energy distributions
        print("Computing Earth Mover's Distance metric on energy...")
        energy_dataset = vectors_from_alignn(dataset, output="energy")
        energy_gen = vectors_from_alignn(generated, output="energy")
        ml_metrics["EMD_energy"] = get_emd(energy_dataset, energy_gen)
        print(f"Energy EMD = {ml_metrics['EMD_energy']}")

    if CONFIG.get("EMD_density"):
        # Compute Earth Mover's Distance on density distributions
        print("Computing Earth Mover's Distance metric on density...")
        densities_dataset = get_densities(dataset)
        densities_generated = get_densities(generated)
        ml_metrics["EMD_density"] = get_emd(
            densities_dataset, densities_generated
        )
        print(f"Density EMD = {ml_metrics['EMD_density']}")

    metrics = {"general": general_metrics, "dft": dft_metrics, "ml": ml_metrics}

    print("Writing output file...")

    with open(args.output, "w", encoding="utf-8") as fp:
        json.dump(metrics, fp, indent=4)

    print(f"Output file successfully written at location '{args.output}'.")
    print(json.dumps(metrics, indent=4))


if __name__ == "__main__":
    main()
