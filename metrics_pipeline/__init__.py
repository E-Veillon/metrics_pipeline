#!/usr/bin/python


import os
from collections.abc import Callable
from typing import Any


# Main paths inside pipeline file tree
ROOT = os.path.dirname(__file__)
"""Absolute path to the metrics pipeline main directory."""

UTILSPATH = os.path.join(ROOT, "utils")
"""Absolute path to the pipeline utils directory containing importable features."""

CONFIGPATH = os.path.join(ROOT, "config")
"""Absolute path to the pipeline config directory containing all YAML config files."""

# Arguments parsing for all executable scripts inside this project
def _parse_input_args(
    cmd_line_func: Callable,
    process_args_func: Callable,
    standalone: bool = True,
    **kwargs
) -> dict[str, Any]:
    """
    Gather and process input arguments, either from command-line or external script.
    
    Parameters:
        cmd_line_func (callable):       Function parsing command line arguments into attributes
                                        of a python object (e.g. argparse.Namespace object).

        process_args_func (callable):   Function asserting and processing input arguments.

        standalone (bool):              Whether parsed script is used directly through
                                        command-line (stand-alone script) or in an external
                                        pipeline script.

        kwargs:                         Input arguments values for external script calls.

    Returns:
        Dict of input arguments names and values.
    """
    if standalone: # Direct use of the script through command line
        print("Command line mode detected.")
        args = cmd_line_func()
        args = args.__dict__
    else: # Indirect use of the script as part of a pipeline script in another file
        print("Pipeline mode detected.")
        args = kwargs

    args: dict[str, Any] = process_args_func(args)

    return args
