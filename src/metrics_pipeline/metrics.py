#!/usr/bin/python
"""
Compute various metrics based on previously done computations.
Computable metrics are:
- Validity,
- Viability,
- Symmetry,
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
from collections.abc import Callable

from pymatgen.core import Structure

from metrics_pipeline.core.utils import (
    parse_input_args, check_type, check_num_value,
    ALL_ELTS_CATEGORIES, get_elts_from_symbol_or_z, get_elts_in_categories,
    filter_by_elements,
    generate_genmat_structures, GenMatStructure, GenMatName
)
from metrics_pipeline.core.genmat_io import (
    check_file_or_dir, check_file_format, load_yaml_as_dict,
    JsonLoader, JsonWriter, PathLike, CONFIGPATH,
    VaspParser, VaspExtractor, ExtractMethod, GenMatFile, StructureFile
)
from metrics_pipeline.core.computations.models import get_crystalnn_fingerprints, vectors_from_alignn
from metrics_pipeline.core.computations.local import get_densities
from metrics_pipeline.core.metrics import (
    StructValidity, Viability, SymmetryClassifier, Unicity, SUN,
    Coverage, EMD, RMSD, FrechetDistance, MetricsData
)


def _get_command_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(
        prog=os.path.basename(__file__), description=__doc__,
        formatter_class=ap.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "generated",
        help="Preprocessed JSON file containing all generated structures."
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
            "Preprocessed JSON file containing a set of reference structures for distribution "
            "comparisons. Necessary for: COV-P, COV-R, FAD, EMD(density), EMD(energy)."
        )
    )
    parser.add_argument(
        "-D", "--database",
        help=(
            "Preprocessed JSON file containing structure data from a reference database "
            "of known structures. Necessary for: Novelty part of the S.U.N. metrics."
        )
    )
    parser.add_argument(
        "-ss", "--sun-summary",
        help=(
            "JSON file containing the summary of phase diagrams instability energies. "
            "The file must be located inside the directory where corresponding VASP run "
            "directories are stored for the runs to be found and properly parsed. "
            "Necessary for: Stability part of S.U.N. metrics."
        )
    )
    parser.add_argument(
        "-rs", "--relax-summary",
        help=(
            "JSON file containing a summary for the relaxation step. "
            "The file must be located inside the directory where corresponding VASP run "
            "directories are stored for the runs to be found and properly parsed. "
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
        "--database-mode", action=ap.BooleanOptionalAction, default=False,
        help=(
            "Only used if Novelty metric computation is enabled. Computes Novelty in two steps, "
            "first initializing the metric with references alone, then comparing each tested "
            "structures one by one against it, instead of initializing and computing everything "
            "at the same time. Can be more stable if '--database' is way bigger than the number "
            "of tested structures. Defaults to %(default)s."
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
    parser.add_argument(
        "--detailed", action=ap.BooleanOptionalAction, default=True,
        help=(
            "Additionally to the metrics values summary file, whether to also output a JSON "
            "file listing individual metrics status for each tested structure, making any "
            "further metrics compounding doable without any extra expensive computations."
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
    args_dict.setdefault("database_mode", False)
    args_dict.setdefault("threshold", 0.4)
    args_dict.setdefault("unique_dataset", True)
    args_dict.setdefault("remove_elts", [])
    args_dict.setdefault("remove_elt_categories", [])
    args_dict.setdefault("detailed", True)

    # Assert set arguments conformity
    check_file_or_dir(args_dict["generated"], "file", allowed_formats="json")
    config_file = os.path.join(CONFIGPATH, args_dict["config"])
    check_file_or_dir(config_file, "file", allowed_formats=("yml","yaml"))

    if args_dict.get("dataset") is not None:
        check_file_or_dir(args_dict["dataset"], "file", allowed_formats="json")

    if args_dict.get("database") is not None:
        check_file_or_dir(args_dict["database"], "file", allowed_formats="json")

    if args_dict.get("sun_summary") is not None:
        check_file_or_dir(args_dict["sun_summary"], "file", allowed_formats="json")

    if args_dict.get("relax_summary") is not None:
        check_file_or_dir(args_dict["relax_summary"], "file", allowed_formats="json")

    check_file_format(args_dict["output"], allowed_formats="json")

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

    check_type(args_dict["database_mode"], "database_mode", (bool,))
    check_type(args_dict["unique_dataset"], "unique_dataset", (bool,))
    check_type(args_dict["detailed"], "detailed", (bool,))

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
    print(f"Stability, Unicity, Novelty (SUN): {config.get('SUN')}")
    print(f"- Stability: {config.get('Stability')}")
    print(f"- Unicity: {config.get('Unicity')}")
    print(f"- Novelty: {config.get('Novelty')}")
    print(f"Avg. Root Mean Square Displacement (RMSD): {config.get('RMSD')}")
    print(f"Coverage - Precision (COV-P): {config.get('COV-P')}")
    print(f"Coverage - Recall (COV-R): {config.get('COV-R')}")
    print(f"Fréchet ALIGNN Distance (FAD): {config.get('FAD')}")
    print(f"Earth Mover's Distance (EMD) on energy: {config.get('EMD_energy')}")
    print(f"Earth Mover's Distance (EMD) on density: {config.get('EMD_density')}")
    print(" ")
    print("------------------------------")
    print(" ")

def _print_elements_removal(args_dict: dict[str, tp.Any]) -> None:
    print("------------------------------")
    print(" ")
    print(
        "REMOVED ELEMENTS: "
        f"{', '.join(args_dict['remove_elts']) if args_dict['remove_elts'] != [] else None}."
    )
    print(
        "REMOVED CATEGORIES: "
        f"{', '.join(args_dict['remove_elt_categories']) if args_dict['remove_elt_categories'] != [] else None}."
    )
    print(
        "ALL REMOVED ELEMENTS (parsed categories): "
        f"{', '.join(args_dict['forbidden_elts'].keys()) if args_dict['forbidden_elts'] != {} else None}."
    )
    print(" ")
    print("------------------------------")

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
        Preprocessed JSON file containing all generated structures.

    config: str | Path
        Configuration file in Yaml format stating which metrics should be computed.
        The file needs to be in 'metrics_pipeline/config' to be found. Defaults to
        metrics_defaults.yaml. This argument adds flexibility by allowing to copy
        the default config file to make several computations presets as needed and
        provide the one needed here.

    dataset: str | Path
        Preprocessed JSON file containing a set of reference structures for distribution
        comparisons. Necessary for: COV-P, COV-R, FAD, EMD(density), EMD(energy).

    database: str | Path
        Preprocessed JSON file containing structure data from a reference database
        of known structures. Necessary for: Novelty part of the S.U.N. metrics.

    sun_summary: str | Path
        JSON file containing the summary of phase diagrams instability energies.
        The file must be located inside the directory where corresponding VASP run
        directories are stored for the runs to be found and properly parsed.
        Necessary for: Stability part of S.U.N. metrics.

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

    database_mode: bool
        Only used if Novelty metric computation is enabled. Computes Novelty in two steps,
        first initializing the metric with references alone, then comparing each tested
        structures one by one against it, instead of initializing and computing everything
        at the same time. Can be more stable if '--database' is way bigger than the number
        of tested structures. Defaults to False.

    remove_elts: list[str], optional
        Structures containing given elements will be removed from reference dataset before
        computing metrics. Saves computation time on metrics if some reference structures
        contain elements that were removed from the generation at preprocessing.
        Elements to track can be specified by their symbol, atomic number, or a mix of both.

    remove_elt_categories: list[str], optional
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

    detailed: bool
        Additionally to the metrics values summary file, whether to also output a JSON
        file listing individual metrics status for each tested structure, making any
        further metrics compounding doable without any extra expensive computations.
    """
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    CONFIG = load_yaml_as_dict(os.path.join(CONFIGPATH, args["config"]), on_error='raise')
    CONFIG.update(CONFIG.pop("SUN", {})) # Add SUN sub-metrics for easier access
    # SUN section activation boolean
    CONFIG["SUN"] = any(CONFIG.get(metric, False) for metric in ("Stability", "Unicity", "Novelty"))
    dataset_needed: bool = (
        CONFIG.get("COV-P", False)
        or CONFIG.get("COV-R", False)
        or CONFIG.get("FAD", False)
        or CONFIG.get("EMD_energy", False)
        or CONFIG.get("EMD_density", False)
    )
    database_needed: bool = CONFIG.get("Novelty", False)
    sun_summary_needed: bool = CONFIG.get("Stability", False)
    relax_summary_needed: bool = CONFIG.get("RMSD", False)

    _print_metrics_config(CONFIG)

    print("===== LOAD NECESSARY DATA FILES =====")

    print("Loading generated structures...")
    gen_file = GenMatFile.from_file(args["generated"])
    generated = gen_file.parse_structures()

    print("Generated structures loaded.")

    dataset = []
    if _match_file_arg_need("dataset", args.get("dataset", None), dataset_needed):
        _print_elements_removal(args)
        print("Loading dataset...")
        dataset = StructureFile.from_file(args["dataset"]).parse_structures()
        dataset, nbr_discarded = filter_by_elements(
            dataset, list(args["forbidden_elts"].keys())
        )
        print(f"{nbr_discarded} reference structures containing forbidden elements were removed.")
        print("Dataset loaded.")

        # remove duplicate structures from the dataset
        if args["unique_dataset"]:
            dataset = generate_genmat_structures(dataset)
            dataset = [
                struct.as_structure() for struct in Unicity(
                    dataset, workers=args.get("workers")
                ).unique_structs
            ]

    database: list[Structure] = []
    if _match_file_arg_need("database", args.get("database", None), database_needed):
        _print_elements_removal(args)
        print("Loading database file...")
        database = StructureFile.from_file(args["database"]).parse_structures()
        database, nbr_discarded = filter_by_elements(
            database, list(args["forbidden_elts"].keys())
        )
        print(f"{nbr_discarded} reference structures containing forbidden elements were removed.")
        print("Database loaded.")

    metastable_names: list[str] = []
    stable_names: list[str] = []
    if _match_file_arg_need("sun-summary", args.get("sun_summary", None), sun_summary_needed):
        print("Loading stability summary file...")
        sun_summary = JsonLoader(args["sun_summary"]).load_as_dict()
        metastable_names = [name for name, s_data in sun_summary.items() if s_data["metastable"]]
        stable_names = [name for name, s_data in sun_summary.items() if s_data["stable"]]
        print("Stability data loaded.")

    relax_structures: list[tuple[GenMatStructure, Structure]] = []
    if _match_file_arg_need("relax-summary", args.get("relax_summary", None), relax_summary_needed):
        print("Loading relaxations summary file...")
        base_dir, summary_name = os.path.split(args["relax_summary"])
        relax_structures_dict = VaspExtractor(
            vasp_parser=VaspParser(base_dir, workers=args.get("workers")),
            method=ExtractMethod.RELAXATION,
            summary_name=summary_name,
            summary_key="converged",
            workers=args.get("workers")
        ).get_data()

        if not relax_structures_dict:
            _warn_summary_location(args["relax_summary"])

        for name, s_data in relax_structures_dict.items():
            init_struct = GenMatStructure(GenMatName(name), s_data["in_struct"])
            final_struct = s_data["out_struct"]
            final_struct.properties["header"] = name
            relax_structures.append((init_struct, final_struct))

        del relax_structures_dict # Free some memory
        print("Data converted.")

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
            "num_metastable", "percent_metastable",
            "num_stable", "percent_stable",
            "num_MSUN", "percent_MSUN",
            "num_SUN", "percent_SUN",
            "RMSD"
        )
    )
    ml_metrics = dict.fromkeys(
        ("precision", "recall", "frechet_distance", "EMD_energy", "EMD_density")
    )

    print("\n===== COMPUTE ACTIVATED METRICS =====")

    # total number of generated structures
    num_generated = len(gen_file)
    round_digits = len(str(num_generated)) - 2
    general_metrics["num_generated"] = num_generated
    general_metrics["num_not_loaded"] = num_generated - len(generated)

    # Generate metrics details dict
    metrics_details: dict[str, MetricsData] | None = {} if args["detailed"] else None
    if metrics_details is not None:
        for struct in generated:
            metrics_details[struct.name] = MetricsData(struct)
    
    def _update_metrics_details(new_data: tp.Iterable[MetricsData]) -> None:
        """Update the detailed dict with newly computed data."""
        if metrics_details is None:
            return
        for data in new_data:
            metrics_details[data.name] |= data


    if CONFIG.get("Validity", False):
        # Validity metric
        print("Computing Validity metric...")
        validity = StructValidity(generated)
        num_valid = len(validity.valid_structs)
        percent_valid = round(num_valid / num_generated * 100, round_digits)
        general_metrics["num_valid"] = num_valid
        general_metrics["percent_valid"] = percent_valid
        print(f"Validity = {percent_valid}%")

        _update_metrics_details(validity.computed_data)

    if CONFIG.get("Viability", False):
        # Viability metric
        print("Computing Viability metric...")
        viability = Viability(generated)
        num_viable = len(viability.viable_structs)
        percent_viable = round(num_viable / num_generated * 100, round_digits)
        general_metrics["num_viable"] = num_viable
        general_metrics["percent_viable"] = percent_viable
        print(f"Viability = {percent_viable}%")

        _update_metrics_details(viability.computed_data)

    if CONFIG.get("Symmetry", False):
        # Symmetry metric
        print("Computing Symmetry metric...")
        symmetrizer = SymmetryClassifier(generated, workers=args.get("workers"))
        num_symmetric = len([data.is_symmetric for data in symmetrizer.computed_data])
        percent_symmetric = round(num_symmetric / num_generated * 100, round_digits)
        general_metrics["num_symmetric"] = num_symmetric
        general_metrics["percent_symmetric"] = percent_symmetric
        print(f"Symmetry = {general_metrics.get('percent_symmetric')}%")

        _update_metrics_details(symmetrizer.computed_data)

    if CONFIG.get("SUN", False):
        # S.U.N. metrics
        print("Computing S.U.N. metrics...")
        # Define all necessary computing flags combinations for conditions readability
        compute_stability = CONFIG.get("Stability", False)
        compute_unicity = CONFIG.get("Unicity", False)
        compute_novelty = CONFIG.get("Novelty", False)
        compute_unique_novel = compute_unicity and compute_novelty
        compute_unmatchable = compute_unicity or compute_novelty
        compute_sun = compute_stability and compute_unicity and compute_novelty

        # We already have the stability summary, thus we deactivate stability computing
        sun_wo_stability = SUN(
            generated, database,
            compute_stability=False,
            compute_unicity=compute_unicity,
            compute_novelty=compute_novelty,
            workers=args.get("workers"),
            database_mode=args["database_mode"]
        )

        if compute_unmatchable:
            num_unmatchable = len(sun_wo_stability.get_computed_subset(unmatchable=True))
            percent_unmatchable = round(num_unmatchable / num_generated * 100, round_digits)
            dft_metrics["num_unmatchable"] = num_unmatchable
            dft_metrics["percent_unmatchable"] = percent_unmatchable

        if compute_unicity:
            num_unique = len(sun_wo_stability.get_computed_subset(unique=True, unmatchable=False))
            percent_unique = round(num_unique / num_generated * 100, round_digits)
            dft_metrics["num_unique"] = num_unique
            dft_metrics["percent_unique"] = percent_unique

        if compute_novelty:
            num_novel = len(sun_wo_stability.get_computed_subset(novel=True, unmatchable=False))
            percent_novel = round(num_novel / num_generated * 100, round_digits)
            dft_metrics["num_novel"] = num_novel
            dft_metrics["percent_novel"] = percent_novel

        if compute_unique_novel:
            unique_novel = sun_wo_stability.get_computed_subset(unique=True, novel=True, unmatchable=False)
            num_unique_novel = len(unique_novel)
            percent_unique_novel = round(num_unique_novel / num_generated * 100, round_digits)
            dft_metrics["num_unique_novel"] = num_unique_novel
            dft_metrics["percent_unique_novel"] = percent_unique_novel

        computed_data = sun_wo_stability.computed_data

        if compute_stability:
            num_metastable = len(metastable_names)
            num_stable = len(stable_names)
            percent_metastable = round(num_metastable / num_generated * 100, round_digits)
            percent_stable = round(num_stable / num_generated * 100, round_digits)
            dft_metrics["num_metastable"] = num_metastable
            dft_metrics["percent_metastable"] = percent_metastable
            dft_metrics["num_stable"] = num_stable
            dft_metrics["percent_stable"] = percent_stable

            for data in computed_data:
                data.is_metastable = data.name in metastable_names
                data.is_stable = data.name in stable_names

        _update_metrics_details(computed_data)

        if compute_sun:
            # Predicates return False if any attribute is None (should never happen)
            msun_predicate: Callable[[MetricsData], bool] = lambda data: (
                bool(data.is_metastable) and bool(data.is_novel) and bool(data.is_unique)
            )
            sun_predicate: Callable[[MetricsData], bool] = lambda data: (
                bool(data.is_stable) and bool(data.is_novel) and bool(data.is_unique)
            )
            num_MSUN = len([data for data in computed_data if msun_predicate(data)])
            num_SUN = len([data for data in computed_data if sun_predicate(data)])
            percent_MSUN = round(num_MSUN / num_generated * 100, round_digits)
            percent_SUN = round(num_SUN / num_generated * 100, round_digits)
            dft_metrics["num_MSUN"] = num_MSUN
            dft_metrics["percent_MSUN"] = percent_MSUN
            dft_metrics["num_SUN"] = num_SUN
            dft_metrics["percent_SUN"] = percent_SUN

        for key, val in dft_metrics.items():
            if key == "RMSD" or val is None:
                continue
            print(f"{key} = {val}")

    if CONFIG.get("RMSD", False):
        # RMSD metric
        print("Computing RMSD metric...")
        _in_structs, _out_structs = zip(*relax_structures)
        in_structs: list[GenMatStructure] = list(_in_structs)
        out_structs: list[Structure] = list(_out_structs)
        rmsd = RMSD(in_structs, out_structs)
        dft_metrics["RMSD"] = round(rmsd.average_rmsd, round_digits)
        print(f"RMSD = {dft_metrics['RMSD']}")

        _update_metrics_details(rmsd.computed_data)

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

        _update_metrics_details(coverage.computed_data)

    if CONFIG.get("FAD", False):
        # Compute Fréchet ALIGNN Distance
        print("Computing Fréchet ALIGNN Distance metric...")
        frechet_distance = FrechetDistance(
            structures=generated,
            ref_structs=dataset,
            transform=vectors_from_alignn,
            device="cpu",
            output="latent"
        )
        ml_metrics["frechet_distance"] = frechet_distance.computed_distance
        print(f"Frechet Distance = {frechet_distance.computed_distance}")

        _update_metrics_details(frechet_distance.computed_data)

    if CONFIG.get("EMD_energy", False):
        # Compute Earth Mover's Distance on energy distributions
        print("Computing Earth Mover's Distance metric on energy...")
        emd_energy = EMD(
            structures=generated,
            ref_structs=dataset,
            transform=vectors_from_alignn,
            device="cpu",
            output="energy"
        )
        ml_metrics["EMD_energy"] = emd_energy.computed_distance
        print(f"Energy EMD = {emd_energy.computed_distance}")

        _update_metrics_details(emd_energy.computed_data)

    if CONFIG.get("EMD_density", False):
        # Compute Earth Mover's Distance on density distributions
        print("Computing Earth Mover's Distance metric on density...")
        emd_density = EMD(
            structures=generated,
            ref_structs=dataset,
            transform=get_densities
        )
        ml_metrics["EMD_density"] = emd_density.computed_distance
        print(f"Density EMD = {emd_density.computed_distance}")

        _update_metrics_details(emd_density.computed_data)

    metrics = {"general": general_metrics, "dft": dft_metrics, "ml": ml_metrics}

    print("Writing output file...")

    JsonWriter(args["output"], metrics, indent=4).write_as_dict()

    print(f"Output file successfully written at location {args['output']}.")
    print(json.dumps(metrics, indent=4))

    if metrics_details is not None:
        print("Writing metrics details file...")
        def _parse_metrics_data(data: MetricsData) -> dict[str, tp.Any]:
            """Remove heavy and useless structure details but keep the name."""
            dct = data.as_dict()
            dct["name"] = dct["structure"]["name"]
            dct.pop("structure", None)
            return dct

        results_dict = {name: _parse_metrics_data(data) for name, data in metrics_details.items()}
        details_output = args["output"].replace(".json", "_details.json")
        JsonWriter(details_output, results_dict, indent=0).write_as_dict()
        print(f"Details file successfully written at location {details_output}.")


if __name__ == "__main__":
    main()
