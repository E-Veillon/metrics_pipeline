#!/usr/bin/python
"""
A script to generate randomly perturbated structures from structures in a CIF file.
"""

import os
import argparse


def assert_args(args: argparse.Namespace) -> None:

    valid_parameters = (
        "sites", "lengths", "angles", 
        "a", "b", "c", "alpha", "beta", "gamma"
    )

    assert os.path.isfile(args.input_file), f"{args.input_file}: no such file found."

    assert args.input_file.endswith(".cif"), "Input file must be of CIF format (.cif)."
    
    if args.output is not None:
        assert args.output.endswith(".cif"), "Output file must be of CIF format (.cif)."

    assert all([param in valid_parameters for param in args.parameters]), (
        "Some of the parameters given are not valid. See --help for a list of supported parameters."
    )
    assert len(args.parameters) == len(set(args.parameters)), (
        "Some parameters are passed more than one time, please verify your --parameters argument."
    )

    if "lengths" in args.parameters:
        assert all([param not in valid_parameters[5:8] for param in args.parameters]), (
            "'lengths' parameter cannot be combined with 'a', 'b' or 'c' parameters."
        )

    if "angles" in args.parameters:
        assert all([param not in valid_parameters[8:] for param in args.parameters]), (
            "'angles' parameter cannot be combined with 'alpha', 'beta' or 'gamma' parameters."
        )

    if args.smallest_perturb is not None:
        assert len(args.smallest_perturb) == len(args.parameters), (
            "The number of values in --smallest-scale argument must match the number of parameters in --parameters."
        )
        assert all([val >= 0.0 for val in args.smallest_perturb]), (
            "Minimal amplitudes must be positive or zero "
            "(whether the perturbation will be positive or negative is random)."
        )

    if args.biggest_perturb is not None:
        assert len(args.biggest_perturb) == len(args.parameters), (
            "The number of values in --biggest-scale argument must match the number of parameters in --parameters."
        )
        assert all([val > 0.0 for val in args.biggest_perturb]), (
            "Maximal amplitudes must be stricly positive "
            "(whether the perturbation will be positive or negative is random)."
        )
    
    assert args.sample_size >= 1, "'sample_size' argument must be strictly positive."

    assert args.workers >= 1, "'workers' argument must be strictly positive."


def main() -> None:
    
    # ARGUMENTS PARSING BLOCK

    parser = argparse.ArgumentParser(
        prog="perturbator.py",
        description=(
            "A script to generate randomly perturbated structures "
            "from given structures in a CIF file."
        ),
        epilog=(
            "Usage example:\n"
            "python perturbator.py structs.cif -p sites a b angles -b 0.05 0.1 0.1 2.0\n"
            "This set will perturb randomly:\n"
            "- sites positions from 0 to 0.05 angstroms from their original positions,\n"
            "- length of lattice vector 'a' from 0 to 0.1 angstroms from its original value,\n"
            "- length of lattice vector 'b' from 0 to 0.1 angstroms from its original value,\n"
            "- value of all lattice angles from 0 to 2.0 degrees from their original values.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the CIF file containing reference structures.\n"
    )
    parser.add_argument(
        "-o","--output",
        type=str,
        default=None,
        help=(
            "Path to the output CIF file to write perturbed structures.\n"
            "By default, it is written at the same path as the input file "
            "with a '_perturb' suffix in name.\n"
        )
    )
    parser.add_argument(
        "-p", "--parameters",
        nargs="*",
        type=str,
        default=("sites",),
        help=(
            "Defines which parameters will be perturbed. Supported args are listed below:\n"
            "- 'sites': perturbs sites positions in lattice (in angstroms).\n"
            "- 'lengths': perturbs all lattice vectors lengths in angstroms "
            "(a, b, c) with same limit amplitudes.\n"
            "- 'angles': perturbs all lattice angles in degrees (alpha, beta, gamma) "
            "with same limit amplitudes.\n"
            "- 'a', 'b', 'c': perturbs corresponding lattice vector length in angstroms.\n"
            "- 'alpha', 'beta', 'gamma': perturbs corresponding lattice angle in degrees.\n"
            "Several args can be given at once in any order, they will be combined to perturb "
            "specific sets of parameters\n"
            "(e.g. '-p sites a b alpha' will perturb sites positions, a, b and alpha "
            "lattice parameters).\n"
            "If some parameters are overlapping, an error is raised.\n"
            "By default, only sites will be perturbed and lattice will remain untouched.\n"
        )
    )
    parser.add_argument(
        "-s", "--smallest-perturb",
        nargs="*",
        type=float,
        default=None,
        help=(
            "Defines the minimal amplitude of applied perturbations\n"
            "(in angstroms for sites positions and lattice vectors, in degrees for lattice angles).\n"
            "Each given value is assigned to parameters in the order given in --parameters.\n"
            "By default, it is set to 0.0 angstrom or 0.0 degree for all given parameters.\n"
        )
    )
    parser.add_argument(
        "-b", "--biggest-perturb",
        nargs="*",
        type=float,
        default=None,
        help=(
            "Defines the maximal amplitude of applied perturbations\n"
            "(in angstroms for sites positions and lattice vectors, in degrees for lattice angles).\n"
            "Each given value is assigned to parameters in the order given in --parameters.\n"
            "By default, it is set to 0.1 angstrom or 1.0 degree for all given parameters.\n"
        )
    )
    parser.add_argument(
        "-n", "--sample-size",
        type=int,
        default=1,
        help="Number of randomly perturbed structures to generate for each input structure.\n"
    )
    parser.add_argument(
        "-r", "--retries",
        default=4,
        help=(
            "Number of times lattice perturbations can be retried "
            "when resulting in an unphysical lattice generation "
            "(i.e. a, b or c < 2*Bohr radius = ~1.06 angstroms, "
            "or volume < (2 * Bohr radius)^3 = ~1.2 angstroms^3) before raising an error."
        )
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to spawn for parallelized steps.\n",
        metavar="int",
    )

    args = parser.parse_args()

    assert_args(args)

    # Manage arguments values and set defaults
    infile = args.input_file

    if args.output is None:
        outfile = infile.replace(".cif", "_perturb.cif")
    else:
        outfile = args.output

    parameters = args.parameters

    if args.smallest_perturb is None:
        min_perturbs = [0.0] * len(parameters)
    else:
        min_perturbs =  args.smallest_perturb

    if args.biggest_perturb is None:
        max_perturbs = []
        for param in parameters:
            match param:
                case "sites"|"lengths"|"a"|"b"|"c":
                    max_perturbs.append(0.1)
                case "angles"|"alpha"|"beta"|"gamma":
                    max_perturbs.append(1.0)
    else:
        max_perturbs = args.biggest_perturb

    sample_size = args.sample_size
    workers = args.workers


    # MAIN BLOCK

    from perturb_utils import (
        get_perturbs_dict, simply_read_cif, write_cif, batch_generate_perturbed_structs
    )

    perturbs_dict = get_perturbs_dict(parameters, min_perturbs, max_perturbs)

    lattice_params = ("a", "b", "c", "alpha", "beta", "gamma")
    modified_lattice = any([perturbs_dict[param] is not None for param in lattice_params])
    
    structures = simply_read_cif(filename=infile, workers=workers)

    perturbed_structs = batch_generate_perturbed_structs(
        structures=structures,
        workers=workers,
        perturbs_dict=perturbs_dict,
        sample_size=sample_size,
        modified_lattice=modified_lattice,
        lattice_retries=args.retries
    )

    write_cif(filename=outfile, structures=perturbed_structs, workers=workers)

if __name__ == "__main__":
    main()