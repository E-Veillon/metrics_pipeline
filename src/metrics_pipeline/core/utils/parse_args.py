"""Parse arguments either from CLI or higher level script for all executable scripts of GenMat."""

import typing as tp
from collections.abc import Callable


def parse_input_args(
    cmd_line_func: Callable,
    process_args_func: Callable,
    standalone: bool = True,
    **kwargs
) -> dict[str, tp.Any]:
    """
    Gather and process input arguments, either from command-line or external script.
    
    Parameters
    ----------
    cmd_line_func: callable
        Function parsing command line arguments into attributes of a python object
        (e.g. argparse.Namespace object).

    process_args_func: callable
        Function asserting and processing input arguments.

    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    kwargs: Any
        Input arguments values for external script calls.

    Returns
    -------
    dict[str, Any]
        Dict of input arguments names and values.
    """
    if standalone: # Direct use of the script through command line
        print("Command line mode detected.")
        args = cmd_line_func()
        args = args.__dict__
    else: # Indirect use of the script as part of a pipeline script in another file
        print("Pipeline mode detected.")
        args = kwargs

    args: dict[str, tp.Any] = process_args_func(args)

    return args