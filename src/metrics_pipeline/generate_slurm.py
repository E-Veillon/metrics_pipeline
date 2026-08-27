#!/usr/bin/python
"""Generate slurm scripts at given paths."""

import os
from pathlib import Path
import argparse as ap
import typing as tp

from core.utils import parse_input_args, check_type, check_num_value
from core.genmat_io import SlurmWriter, check_file_or_dir


def _get_cmd_line_args() -> ap.Namespace:
    parser = ap.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "--output", "-o",
        help="Path to write the slurm script to."
    )
    parser.add_argument(
        "--filename", "-f", default="run.slurm",
        help="Optional name to give to the slurm script. Defaults to '%(default)s'."
    )
    parser.add_argument(
        "--jobname", default="jobname",
        help="Name of the job in the job queue. Defaults to '%(default)s'."
    )
    parser.add_argument(
        "--slurm-output",
        help=(
            "Name of the slurm output file capturing stdout output. "
            "Defaults to 'slurm_%j.out', or 'slurm_%A_%a.out' if 'array' is given. "
            "See Slurm documentation for more details about '%' shortcuts behaviors."
        )
    )
    parser.add_argument(
        "--slurm-error",
        help=(
            "Name of the slurm error file capturing stderr output."
            "Defaults to 'slurm_%j.err', or 'slurm_%A_%a.err' if 'array' is given. "
            "See Slurm documentation for more details about '%' shortcuts behaviors."
        )
    )
    parser.add_argument(
        "--nodes", type=int, default=1,
        help="Number of reserved computation nodes. Must be at least 1. Defaults to %(default)s."
    )
    parser.add_argument(
        "--ntasks", type=int, default=0,
        help=(
            "Max number of tasks to do in parallel in total. If not given, `ntasks_per_node` "
            "must be given."
        )
    )
    parser.add_argument(
        "--ntasks_per_node", type=int, default=0,
        help=(
            "Max number of tasks to do in parallel in each node. If not given, `ntasks` "
            "must be given."
        )
    )
    parser.add_argument(
        "--n_gpus", type=int, default=0,
        help="Number of reserved GPUs in each node. Defaults to %(default)s for CPU only."
    )
    parser.add_argument(
        "--cpus_per_task", type=int, default=1,
        help="Number of reserved CPU cores for each task. Defaults to %(default)s."
    )
    parser.add_argument(
        "--time",
        help=(
            "Max timeout for the job before stopping. Format is 'd-hh:mm:ss', "
            "where d stands for days, h for hours, m for minutes and s for seconds."
            "If the job is less than a day, the 'd-' part is optional."
        )
    )
    parser.add_argument(
        "--constraint",
        help=(
            "Apply a constraint to the job parameters, such as reserving only specific nodes."
            "Must be a valid constraint defined for the slurm installation. See cluster "
            "documentation for more details."
        )
    )
    parser.add_argument(
        "--partition",
        help="Name of a specific partition to target for the run."
    )
    parser.add_argument(
        "--hint",
        help=(
            "Apply a specific property to the run, such as deactivationg hyperthreading. "
            "See cluster documentation for more details."
        )
    )
    parser.add_argument(
        "--qos",
        help="Apply a defined Quality of Service. See cluster documentation for more details."
    )
    parser.add_argument(
        "--account",
        help="Define an account to debit computation hours from, if applicable."
    )
    parser.add_argument(
        "--array",
        help=(
            "Define a queue of jobs using the same script with indices. "
            "Use the slurm environment variable `$SLURM_ARRAY_TASK_ID` to use the array index "
            "of each queued job and distinguish between each job similar actions."
        )
    )
    parser.add_argument(
        "--content-file",
        help=(
            "Path to a file containing the bash code to copy into generated script."
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
    args_dict.setdefault("output", os.getcwd())
    args_dict.setdefault("filename", "run.slurm")
    args_dict.setdefault("jobname", "jobname")
    args_dict.setdefault("nodes", 1)
    args_dict.setdefault("ntasks", 0)
    args_dict.setdefault("ntaks_per_node", 0)
    args_dict.setdefault("n_gpus", 0)
    args_dict.setdefault("cpus_per_task", 1)

    # Assert set arguments conformity
    check_type(args_dict["output"], "output", (str, Path))
    check_type(args_dict["filename"], "filename", (str,))
    check_type(args_dict["jobname"], "jobname", (str,))

    for arg_name in (
        "slurm_output", "slurm_error", "time", "constraint", "partition",
        "hint", "qos", "account", "array"
    ):
        check_type(args_dict.get(arg_name), arg_name, (str, type(None)))

    for arg_name in ("nodes", "ntasks", "ntasks_per_node", "n_gpus", "cpus_per_task"):
        check_type(args_dict[arg_name], arg_name, (int,))
        check_num_value(args_dict[arg_name], arg_name, ">=", 0)

    assert args_dict["ntasks"] or args_dict["ntasks_per_node"], ValueError(
        "At least one of either 'ntasks' or 'ntasks_per_node' must be > 0."
    )

    # Additional arguments processing
    args_dict["outfile"] = os.path.join(args_dict.pop("output"), args_dict.pop("filename"))

    if args_dict.get("slurm_output") is not None:
        args_dict["output"] = args_dict.pop("slurm_output")

    if args_dict.get("slurm_error") is not None:
        args_dict["error"] = args_dict.pop("slurm_error")

    if args_dict.get("content_file") is not None:
        check_file_or_dir(args_dict["content_file"], "file")
        with open(args_dict["content_file"], "rt", encoding="utf-8") as fp:
            args_dict["content_lines"] = fp.read().splitlines()

        args_dict.pop("content_file")

    return args_dict


def main(standalone: bool = True, **kwargs) -> None:
    """
    Generate slurm scripts at given paths.
    
    Parameters
    ----------
    standalone (bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    jobname: str, optional
        Name of the job in the job queue. Defaults to "jobname".

    slurm_output: str, optional
        Name of the slurm output file capturing stdout output.
        Defaults to 'slurm_%j.out', or 'slurm_%A_%a.out' if 'array' is given.
        See Slurm documentation for more details about '%' shortcuts behaviors.

    slurm_error: str, optional
        Name of the slurm error file capturing stderr output.
        Defaults to 'slurm_%j.err', or 'slurm_%A_%a.err' if 'array' is given.
        See Slurm documentation for more details about '%' shortcuts behaviors.

    nodes: int, optional
        Number of reserved computation nodes. Must be at least 1.
        Defaults to 1.

    ntasks: int, optional
        Max number of tasks to do in parallel in total. If not given,
        `ntasks_per_node` must be given.

    ntasks_per_node: int, optional
        Max number of tasks to do in parallel in each node. If not given,
        `ntasks` must be given.

    n_gpus: int, optional
        Number of reserved GPUs in each node. Defaults to 0 for CPU only execution.

    cpus_per_task: int, optional
        Number of reserved CPU cores for each task. Defaults to 1 for sequential execution.
        
    time: str, optional
        Max timeout for the job before stopping. Format is "d-hh:mm:ss",
        where d stands for days, h for hours, m for minutes and s for seconds.
        If the job is less than a day, the "d-" part is optional.

    constraint: str, optional
        Apply a constraint to the job parameters, such as reserving only specific nodes.
        Must be a valid constraint defined for the slurm installation. See cluster documentation
        for more details.

    partition: str, optional
        Name of a specific partition to target for the run.

    hint: str, optional
        Apply a specific property to the run, such as deactivationg hyperthreading.
        See cluster documentation for more details.

    qos: str, optional
        Apply a defined Quality of Service. See cluster documentation for more details.

    account: str, optional
        Define an account to debit computation hours from, if applicable.

    array: str
        Define a queue of jobs using the same script with indices.
        Use the slurm environment variable `$SLURM_ARRAY_TASK_ID` to use the array index
        of each queued job and distinguish between each job similar actions.

    content_file: str | Path, optional
        "Path to a file containing the bash code to copy into generated script."
    """
    args = parse_input_args(_get_cmd_line_args, _process_input_args, standalone, **kwargs)
    outfile = args.pop("outfile")
    writer = SlurmWriter(**args)
    writer.write_script(outfile)


if __name__ == "__main__":
    main()