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
        help="Threshold for the computation of coverage recall and coverage precision metrics."
    )

    args = parser.parse_args()

    import os
    import json
    from copy import deepcopy

    import numpy as np
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.analysis.structure_matcher import StructureMatcher

    from screening_pipeline.utils.vasp_io import batch_extract_vasp_structures
    from screening_pipeline.utils.cif_io import read_cif
    from screening_pipeline.utils.matcher import batch_group_by_equivalence, remove_equivalent
    from screening_pipeline.utils.ml_vectors import vectors_from_alignn
    from screening_pipeline.utils.distribution import (
        recall,
        precision,
        frechet_distance,
        wasserstein_distance,
    )
    from screening_pipeline.utils.density import get_densities
    from screening_pipeline.utils.rmsd import rmsd_from_structures

    full_generated, _, _ = read_cif(
        filename=args.generated,
        workers=args.workers,
        keep_rare_gases=True,
        keep_rare_earths=True
    )
    if not args.no_rare_gas_check or not args.no_rare_earth_check:
        pruned_generated, _, _ = read_cif(
            filename=args.generated,
            workers=args.workers,
            keep_rare_gases=args.no_rare_gas_check,
            keep_rare_earths=args.no_rare_earth_check
        )
    else:
        pruned_generated = deepcopy(full_generated)

    symmetrized, _, _ = read_cif(
        filename=args.symmetrized,
        workers=args.workers,
        keep_rare_gases=True,
        keep_rare_earths=True
    )
    dataset, _, _ = read_cif(
        filename=args.dataset,
        workers=args.workers,
        keep_rare_gases=True,
        keep_rare_earths=True
    )

    with open(args.summary, "r") as fp:
        summary = json.load(fp)

    vasp_structures = batch_extract_vasp_structures(
        calc_dirs=[struct["path"] for struct in summary],
        workers=args.workers
    )

    # remove duplicate structures from the dataset
    dataset, _ = remove_equivalent(
        structures=dataset, workers=args.workers, keep_equivalent=False
    )

    # total
    num_generated = len(full_generated)
    num_generated_wo_rare = len(pruned_generated)

    # Unique count
    num_unique = len(symmetrized)

    """
    # novel count
    concat_novel, _ = remove_equivalent(
        structures=pruned_generated + dataset, workers=args.workers, keep_equivalent=False
    )
    num_novel = len(concat_novel) - len(dataset)

    # novel + unique count
    concat_novel_unique, _ = remove_equivalent(
        structures=symmetrized + dataset, workers=args.workers, keep_equivalent=False
    )
    num_novel_unique = len(concat_novel_unique) - len(dataset)

    # novel + unique + stable count
    num_novel_unique_stable = sum(map(lambda x: x["stable"], summary))
    # TODO: cette ligne ne garantit pas la nouveauté
    # TODO: Modifier pour obtenir l'intersection entre les structures nouvelles et le résumé
    """

    # Stable + Unique count
    paths_stable_unique = list(filter(lambda data: data["stable"], summary))
    num_stable_unique = len(paths_stable_unique)

    # Stable + Unique + Novel count (SUN)
    structs_stable_unique = list(map(
        lambda data: Poscar.from_file(os.path.join(data["path"], "POSCAR")),
        paths_stable_unique
    ))
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
    num_stable_unique_novel = len(sun_structs)

    # RMSD
    in_structs = [s for s, _ in vasp_structures]
    out_structs = [s for _, s in vasp_structures]
    rmsd = np.mean(rmsd_from_structures(in_structs, out_structs))

    # machine learning
    latent_dataset = vectors_from_alignn(dataset,output="latent")
    latent_gen = vectors_from_alignn(full_generated,output="latent")

    p = precision(latent_gen, latent_dataset, args.threshold)
    r = recall(latent_gen, latent_dataset, args.threshold)
    fd = frechet_distance(latent_gen, latent_dataset)

    energy_dataset = vectors_from_alignn(dataset,output="energy")
    energy_gen = vectors_from_alignn(full_generated,output="energy")
    emd_energy = wasserstein_distance(energy_dataset,energy_gen)

    densities_dataset=get_densities(dataset)
    densities_generated=get_densities(full_generated)
    emd_density = wasserstein_distance(densities_dataset,densities_generated)


    metrics = {
        "dft": {
            "num_generated": num_generated,
            "num_generated_wo_rare":  num_generated_wo_rare,
            "num_unique": num_unique,
            #"num_novel": num_novel,
            #"num_novel_unique": num_novel_unique,
            "num_stable_unique": num_stable_unique,
            "num_stable_unique_novel": num_stable_unique_novel,
            "SUN": num_stable_unique_novel / num_generated_wo_rare,
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
