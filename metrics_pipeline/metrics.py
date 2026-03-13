#!/usr/bin/python
"""
Compute various metrics based on previously done computations.
Computable metrics are:
- Validity,
- Viability,
- Symmetry,
- Elementary Metastability,
- Stability, Unicity, Novelty (S.U.N),
- Average Root Mean Square Displacement (RMSD),
- Coverage - Precision (COV-P) and Coverage - Recall (COV-R),
- Fréchet ALIGNN Distance (FAD),
- Earth Mover's Distance (EMD) on energy and density distributions.
"""

import os
import json
import argparse as ap
import typing as tp
import warnings

from src.utils import (
    parse_input_args, check_type, check_num_value,
    discard_rare_gas_structures, discard_rare_earth_structures,
    ALL_ELTS_CATEGORIES, get_elts_from_symbol_or_z, get_elts_in_categories,
    filter_by_elements
)
from src.io import (
    check_file_or_dir, check_file_format, load_yaml_as_dict,
    JsonWriter, CIFFile, PathLike, CONFIGPATH, VaspParser, VaspExtractor, ExtractMethod
)
from src.computations.models import get_crystalnn_fingerprints, vectors_from_alignn
from src.computations.local import get_densities
from src.metrics import (
    StructValidity, Viability, Symmetry,
    ElementaryMetastability, Novelty, Unicity, SUN,
    Coverage, EMD, RMSD, FrechetDistance
)


def _get_command_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(
        prog=os.path.basename(__file__), description=__doc__,
        formatter_class=ap.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "generated",
        help="Preprocessed CIF file containing all generated structures."
    )
    parser.add_argument(
        "-c", "--config", default="metrics_defaults.yaml",
        help=(
            "Configuration file in Yaml format stating which metrics should be computed. "
            f"The file needs to be in {CONFIGPATH} to be found. "
            "Defaults to %(default)s. This argument adds flexibility by allowing to copy "
            "the default config file to make several computations presets as needed and "
            "provide the one needed here."
        )
    )
    parser.add_argument(
        "-d", "--dataset",
        help=(
            "CIF file containing the list of known structures. "
            "Necessary for: Elementary Metastability, S.U.N., COV-P, COV-R, FAD, EMD(density), EMD(energy)."
        )
    )
    parser.add_argument(
        "-ss", "--sun-summary",
        help=(
            "JSON file containing the summary of phase diagrams instability energies. "
            "Necessary for: S.U.N."
        )
    )
    parser.add_argument(
        "-rs", "--relax-summary",
        help=(
            "JSON file containing a summary for the relaxation step. "
            "Necessary for: RMSD."
        )
    )
    parser.add_argument(
        "-o", "--output", default="metrics.json",
        help="Output file containing the calculated metrics (json format).",
    )
    parser.add_argument(
        "--workers", "-w", type=int, metavar="int",
        help=(
            "Number of processes to use in parallel. If not given, will use default of "
            "`tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()` "
            "and execute sequentially."
        )
    )
    parser.add_argument(
        "-r", "--remove-elts", nargs="*",
        help=(
            "Structures containing given elements will be removed from reference dataset before "
            "computing metrics. Saves computation time on metrics if some reference structures "
            "contain elements that were removed from the generation at preprocessing. "
            "Elements to track can be specified by their symbol, atomic number, or a mix of both."
        )
    )
    parser.add_argument(
        "-R", "--remove-elt-categories", nargs="*",
        help=(
            "Pass valid element categories to eliminate reference dataset structures containing "
            "any element from these categories. Supported categories are: "
            f"{', '.join(sorted(ALL_ELTS_CATEGORIES))}."
        )
    )
    parser.add_argument(
        "-t", "--threshold", default=0.4, type=float,
        help="Threshold for the computation of coverage recall and coverage precision metrics.",
    )
    parser.add_argument(
        "--unique-dataset", action=ap.BooleanOptionalAction, default=True,
        help=(
            "Whether to use StructureMatcher to remove duplicate structures from reference "
            "dataset before computing metrics. Can save computation time on metrics if a lot "
            "of duplicates can be removed, but can also take significant time to perform "
            "matching on a big dataset. Activated by default. Deactivate if your dataset is "
            "unlikely to have a lot of duplicates."
        )
    )
    args = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}

    args_dict.setdefault("config", "metrics_defaults.yaml")
    args_dict.setdefault("output", "metrics.json")
    args_dict.setdefault("threshold", 0.4)
    args_dict.setdefault("unique_dataset", True)
    args_dict.setdefault("remove_elts", [])
    args_dict.setdefault("remove_elt_categories", [])

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("generated"), "file", allowed_formats="cif")
    config_file = os.path.join(CONFIGPATH, args_dict["config"])
    check_file_or_dir(config_file, "file", allowed_formats=("yml","yaml"))

    if args_dict.get("dataset") is not None:
        check_file_or_dir(args_dict.get("dataset"), "file", allowed_formats="cif")

    if args_dict.get("sun_summary") is not None:
        check_file_or_dir(args_dict.get("sun_summary"), "file", allowed_formats="json")

    if args_dict.get("relax_summary") is not None:
        check_file_or_dir(args_dict.get("relax_summary"), "file", allowed_formats="json")

    check_file_format(args_dict.get("output"), allowed_formats="json")

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">=", 0)

    check_type(args_dict["remove_elts"], "remove_elts", (list,))
    for idx, elt in enumerate(args_dict["remove_elts"]):
        check_type(elt, f"remove_elts[{idx}]", (str,))

    check_type(args_dict["remove_elt_categories"], "remove_elt_categories", (list,))
    for idx, elt in enumerate(args_dict["remove_elt_categories"]):
        check_type(elt, f"remove_elt_categories[{idx}]", (str,))
        if not elt in ALL_ELTS_CATEGORIES:
            raise ValueError(f"The category {elt!r} is not supported.")

    check_type(args_dict["threshold"], "threshold", (float,))
    check_num_value(args_dict["threshold"], "threshold", ">", 0.0)

    check_type(args_dict["unique_dataset"], "unique_dataset", (bool,))

    # Parse forbidden elements
    forbidden_elts = get_elts_in_categories(args_dict["remove_elt_categories"])
    forbidden_elts.update(get_elts_from_symbol_or_z(args_dict["remove_elts"]))
    args_dict["forbidden_elts"] = forbidden_elts

    return args_dict


def _match_file_arg_need(
    arg_name: str, filename: tp.Optional[PathLike] = None, is_needed: bool = False
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
    print(f"Viability: {config.get('Viability')}")
    print(f"Symmetry: {config.get('Symmetry')}")
    print(f"Elementary Metastability: {config.get('ElementaryMetastability')}")
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


def _warn_summary_location(summary_file: PathLike) -> None:
    warnings.warn(
        f"No structure data could be extracted from {summary_file}. "
        "Either none of the structures passed the previous filter or the "
        "summary file is not located in the right directory. Make sure "
        "the summary file is located in the directory containing corresponding "
        "VASP run directories."
    )


def main(standalone: bool = True, **kwargs) -> None:
    """
    Compute various metrics based on previously done computations.
    Computable metrics are:
    - Validity,
    - Viability,
    - Symmetry,
    - Elementary Metastability,
    - Stability, Unicity, Novelty (S.U.N),
    - Average Root Mean Square Displacement (RMSD),
    - Coverage - Precision (COV-P) and Coverage - Recall (COV-R),
    - Fréchet ALIGNN Distance (FAD),
    - Earth Mover's Distance (EMD) on energy and density distributions.

    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    generated: str | Path
        Preprocessed CIF file containing all generated structures.

    config: str | Path
        Configuration file in Yaml format stating which metrics should be computed.
        The file needs to be in 'metrics_pipeline/config' to be found. Defaults to
        metrics_defaults.yaml. This argument adds flexibility by allowing to copy
        the default config file to make several computations presets as needed and
        provide the one needed here.

    dataset: str | Path
        Cif file containing the list of known structures. Necessary for: Elementary Metastability, S.U.N.,
        COV-P, COV-R, FAD, EMD(density), EMD(energy).

    sun_summary: str | Path
        JSON file containing the summary of phase diagrams instability energies.
        The file must be located inside the directory where corresponding VASP run
        directories are stored for the runs to be found and properly parsed.
        Necessary for: S.U.N.

    relax_summary: str | Path
        JSON file containing a summary for the relaxation step.
        The file must be located inside the directory where corresponding VASP run
        directories are stored for the runs to be found and properly parsed.
        Necessary for: RMSD.

    output: str | Path
        Output file containing the calculated metrics (json format).

    workers: int
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    remove_elts: list[str], optional
        Structures containing given elements will be removed from reference dataset before
        computing metrics. Saves computation time on metrics if some reference structures
        contain elements that were removed from the generation at preprocessing.
        Elements to track can be specified by their symbol, atomic number, or a mix of both.

    remve_elt_categories: list[str], optional
        Pass valid element categories to eliminate reference dataset structures containing
        any element from these categories.

    threshold: float
        Threshold for the computation of coverage recall and coverage precision metrics.

    unique_dataset: bool
        Whether to use StructureMatcher to remove duplicate structures from reference
        dataset before computing metrics. Can save computation time on metrics if a lot
        of duplicates can be removed, but can also take significant time to perform
        matching on a big dataset. Activated by default. Deactivate if your dataset is
        unlikely to have a lot of duplicates.
    """
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    CONFIG = load_yaml_as_dict(os.path.join(CONFIGPATH, args["config"]), on_error='raise')
    dataset_needed: bool = (
        CONFIG.get("ElementaryMetastability")
        or CONFIG.get("SUN", False)
        or CONFIG.get("COV-P", False)
        or CONFIG.get("COV-R", False)
        or CONFIG.get("FAD", False)
        or CONFIG.get("EMD_energy", False)
        or CONFIG.get("EMD_density", False)
    )
    sun_summary_needed = CONFIG.get("ElementaryMetastability", False) or CONFIG.get("SUN", False)
    relax_summary_needed = CONFIG.get("RMSD", False)

    _print_metrics_config(CONFIG)

    print("===== LOAD NECESSARY DATA FILES =====")

    print("Loading generated structures...")
    gen_file = CIFFile.from_file(
        args["generated"],
        workers=args.get("workers")
    )
    gen_cifs = gen_file.get_cifs()
    generated, _ = gen_file.parse_structures()

    print("Generated structures loaded.")

    if _match_file_arg_need("dataset", args.get("dataset", None), dataset_needed):
        print("Loading dataset...")
        data_file = CIFFile.from_file(
            args["dataset"],
            workers=args.get("workers")
        )
        data_cifs = data_file.get_cifs()
        data_cifs, nbr_discarded = filter_by_elements(
            data_cifs, list(args["forbidden_elts"].keys()), format="cif"
        )
        print(f"{nbr_discarded} reference structures containing forbidden elements were removed.")
        data_file.clear()
        data_file.add_cifs(data_cifs)
        dataset, _ = data_file.parse_structures()

        print("Dataset loaded.")
        # remove duplicate structures from the dataset
        if args["unique_dataset"]:
            dataset = Unicity(dataset, workers=args.get("workers")).unique_structs
    else:
        dataset = []

    if _match_file_arg_need("sun-summary", args.get("sun_summary", None), sun_summary_needed):
        print("Loading stability summary file...")
        base_dir, summary_name = os.path.split(args["sun_summary"])
        stable_structs_pairs = VaspExtractor(
            vasp_parser=VaspParser(base_dir),
            method=ExtractMethod.RELAXATION,
            summary_name=summary_name,
            summary_key="stable",
            workers=args.get("workers")
        ).get_data()

        if not stable_structs_pairs:
            _warn_summary_location(args["sun_summary"])

        stable_structs = [
            (name, s_data["in_struct"]) for name, s_data in stable_structs_pairs.items()
        ]
        del stable_structs_pairs # Free some memory
        print("Data converted.")
    else:
        stable_structs = []

    if _match_file_arg_need("relax-summary", args.get("relax_summary", None), relax_summary_needed):
        print("Loading relaxations summary file...")
        base_dir, summary_name = os.path.split(args["relax_summary"])
        relax_structures_dict = VaspExtractor(
            vasp_parser=VaspParser(base_dir),
            method=ExtractMethod.RELAXATION,
            summary_name=summary_name,
            summary_key="converged",
            workers=args.get("workers")
        ).get_data()

        if not relax_structures_dict:
            _warn_summary_location(args["relax_summary"])

        relax_structures = list(relax_structures_dict.items())
        del relax_structures_dict # Free some memory
        print("Data converted.")
    else:
        relax_structures = []

    general_metrics = dict.fromkeys(
        (
            "num_generated", "num_not_loaded",
            "num_valid", "percent_valid",
            "num_viable", "percent_viable",
            "num_symmetric", "percent_symmetric"
        )
    )
    dft_metrics = dict.fromkeys(
        (
            "num_elem_metastable", "percent_elem_metastable",
            "num_unmatchable", "percent_unmatchable",
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

    print("\n===== COMPUTE ACTIVATED METRICS =====")

    # total number of generated structures
    num_generated = len(gen_cifs)
    round_digits = len(str(num_generated)) - 2
    general_metrics["num_generated"] = num_generated
    general_metrics["num_not_loaded"] = num_generated - len(generated)

    if CONFIG.get("Validity", False):
        # Validity metric
        print("Computing Validity metric...")
        num_valid = len(StructValidity(generated).valid_structs)
        percent_valid = round(num_valid / num_generated * 100, round_digits)
        general_metrics["num_valid"] = num_valid
        general_metrics["percent_valid"] = percent_valid
        print(f"Validity = {percent_valid}%")

    if CONFIG.get("Viability", False):
        # Viability metric
        print("Computing Viability metric...")
        num_viable = len(Viability(generated).viable_structs)
        percent_viable = round(num_viable / num_generated * 100, round_digits)
        general_metrics["num_viable"] = num_viable
        general_metrics["percent_viable"] = percent_viable
        print(f"Viability = {percent_viable}%")

    if CONFIG.get("Symmetry", False):
        # Symmetry metric
        print("Computing Symmetry metric...")
        num_symmetric = len(
            Symmetry(generated, symprec=0.1, workers=args.get("workers")).symmetric_structs
        )
        percent_symmetric = round(num_symmetric / num_generated * 100, round_digits)
        general_metrics["num_symmetric"] = num_symmetric
        general_metrics["percent_symmetric"] = percent_symmetric
        print(f"Symmetry = {general_metrics.get('percent_symmetric')}")

    # FIXME: get actual energies to compute this, unavailable until then
    if CONFIG.get("ElementaryMetastability", False) and False:
        ref_unaries = [struct for struct in dataset if struct.composition.is_element]
        metastability = ElementaryMetastability(generated, ref_unaries)
        dft_metrics["num_elem_metastable"] = len(metastability.metastable_structs)
        prop_metastables = len(metastability.metastable_structs) / num_generated
        dft_metrics["percent_elem_metastable"] = round(prop_metastables * 100, 6)
        print(f"Elem. Metastability = {dft_metrics.get('percent_elem_metastable')}")

    if CONFIG.get("SUN", False):
        # S.U.N. metrics
        print("Computing S.U.N. metrics...")
        sun = SUN(generated, dataset, workers=args.get("workers"))

        num_unmatchable = len(sun.get_computed_subset(unmatchable=True))
        percent_unmatchable = round(num_unmatchable / num_generated * 100, round_digits)
        num_unique = len(sun.get_computed_subset(unique=True, unmatchable=False))
        percent_unique = round(num_unique / num_generated * 100, round_digits)
        num_novel = len(sun.get_computed_subset(novel=True, unmatchable=False))
        percent_novel = round(num_novel / num_generated * 100, round_digits)
        num_unique_novel = len(sun.get_computed_subset(unique=True, novel=True, unmatchable=False))
        percent_unique_novel = round(num_unique_novel / num_generated * 100, round_digits)
        num_stable = len(sun.get_computed_subset(stable=True))
        percent_stable = round(num_stable / num_generated * 100, round_digits)
        num_SUN = len(sun.get_computed_subset(stable=True, unique=True, novel=True, unmatchable=False))
        percent_SUN = round(num_SUN / num_generated * 100, round_digits)

        dft_metrics["num_unmatchable"] = num_unmatchable
        dft_metrics["percent_unmatchable"] = percent_unmatchable
        dft_metrics["num_unique"] = num_unique
        dft_metrics["percent_unique"] = percent_unique
        dft_metrics["num_novel"] = num_novel
        dft_metrics["percent_novel"] = percent_novel
        dft_metrics["num_unique_novel"] = num_unique_novel
        dft_metrics["percent_unique_novel"] = percent_unique_novel
        dft_metrics["num_stable"] = num_stable
        dft_metrics["percent_stable"] = percent_stable
        dft_metrics["num_SUN"] = num_SUN
        dft_metrics["percent_SUN"] = percent_SUN

        for key, val in dft_metrics.items():
            if key == "RMSD":
                continue
            print(f"{key} = {val}")

    if CONFIG.get("RMSD", False):
        # RMSD metric
        print("Computing RMSD metric...")
        in_structs = [s_data["in_struct"] for _, s_data in relax_structures]
        out_structs = [s_data["out_struct"] for _, s_data in relax_structures]
        dft_metrics["RMSD"] = round(RMSD(in_structs, out_structs).average_rmsd, round_digits)
        print(f"RMSD = {dft_metrics['RMSD']}")

    if CONFIG.get("COV-P", False) or CONFIG.get("COV-R", False):
        # Compute Coverage (Precision, Recall)
        print("Computing fingerprints for Coverage metrics...")
        coverage = Coverage(
            generated, dataset,
            transform=get_crystalnn_fingerprints,
            compute_precision=CONFIG.get("COV-P", False),
            compute_recall=CONFIG.get("COV-R", False),
            threshold=args["threshold"],
            workers=args.get("workers")
        )
        ml_metrics["precision"] = coverage.precision
        ml_metrics["recall"] = coverage.recall

    if CONFIG.get("FAD", False):
        # Compute Fréchet ALIGNN Distance
        print("Computing Fréchet ALIGNN Distance metric...")
        frechet_distance = FrechetDistance(
            structures=generated,
            ref_structs=dataset,
            transform=vectors_from_alignn,
            output="latent"
        ).computed_distance
        ml_metrics["frechet_distance"] = frechet_distance
        print(f"Frechet Distance = {frechet_distance}")

    if CONFIG.get("EMD_energy", False):
        # Compute Earth Mover's Distance on energy distributions
        print("Computing Earth Mover's Distance metric on energy...")
        emd_energy = EMD(
            structures=generated,
            ref_structs=dataset,
            transform=vectors_from_alignn,
            output="energy"
        ).computed_distance
        ml_metrics["EMD_energy"] = emd_energy
        print(f"Energy EMD = {emd_energy}")

    if CONFIG.get("EMD_density", False):
        # Compute Earth Mover's Distance on density distributions
        print("Computing Earth Mover's Distance metric on density...")
        emd_density = EMD(
            structures=generated,
            ref_structs=dataset,
            transform=get_densities
        ).computed_distance
        ml_metrics["EMD_density"] = emd_density
        print(f"Density EMD = {emd_density}")

    metrics = {"general": general_metrics, "dft": dft_metrics, "ml": ml_metrics}

    print("Writing output file...")

    JsonWriter(args["output"], metrics, indent=4).write_as_dict()

    print(f"Output file successfully written at location {args['output']}.")
    print(json.dumps(metrics, indent=4))


if __name__ == "__main__":
    main()
