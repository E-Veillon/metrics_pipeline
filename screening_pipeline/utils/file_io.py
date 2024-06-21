#!/usr/bin/python
"""
Implements functions to manage and operate on paths.
"""


import os
import warnings
from typing import Literal
from pathlib import Path
from ruamel.yaml import YAML

# LOCAL IMPORTS
from custom_types import PathLike

# Main paths inside pipeline file tree
MAINDIRPATH = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
"""Absolute path to the main pipeline directory."""
SCRIPTSPATH = os.path.join(MAINDIRPATH, "scripts")
"""Absolute path to the pipeline scripts directory containing executable scripts."""
UTILSPATH   = os.path.join(MAINDIRPATH, "utils")
"""Absolute path to the pipeline utils directory containing importable features."""
CONFIGPATH  = os.path.join(MAINDIRPATH, "config")
"""Absolute path to the pipeline config directory containing all VASP config YAML files."""


########################################


def check_file_or_dir(
    path: PathLike,
    file_or_dir: Literal["file", "dir"] = "file",
    *,
    format: str|None = None
) -> None:
    """
    Verify existence and optionally extension format of given path.
    
    Parameters:
        path (str|Path):    Path to verify.

        file_or_dir (str):  Whether the path should lead to a file or a directory.
                            If the path exists but is not the right data type,
                            an error will still be raised for not finding it.

        format (str):       If the path should lead to a file with a specific format
                            extension, provide here wanted extension without the dot
                            separator (e.g. "txt" and not ".txt").
    """
    assert isinstance(path, (str, Path))
    assert file_or_dir in {"file", "dir"}
    assert isinstance(format, str) or format is None

    path = str(path)

    if file_or_dir == "dir" and not os.path.isdir(path):
        raise FileNotFoundError(
            f"{path}: No such directory found."
        )
    if file_or_dir == "file" and not os.path.isfile(path):
        raise FileNotFoundError(
            f"{path}: No such file found."
        )
    if (
        file_or_dir == "file"
        and format is not None
        and not path.endswith("." + format)
    ):
        file_ext = path.split(sep=".")[-1]
        raise ValueError(
            f"{path}: expected file format is '{format}', "
            f"got '{file_ext}' format instead."
        )


########################################


def add_new_dir(base_dir: PathLike, *new_dirs: str) -> str:
    """
    Creates a new path of sub-directories inside given base directory.
    If part of the path already exists, only lacking subdirs are created.
    The behaviour is like os.makedirs, but the complete path is returned
    as a string once created.
    
    Parameters:
        base_dir (str|Path):    The base directory inside which
                                the new one will be created.

        *dirs (str):            The name of the new subdirectory to create.

    Returns:
        str: path pointing to the new subdirectory.
    """

    assert isinstance(base_dir, (str, Path))
    assert os.path.isdir(str(base_dir))
    assert all(isinstance(dir, str) for dir in new_dirs)

    new_path = os.path.join(str(base_dir), *new_dirs)
    os.makedirs(new_path, exist_ok=True)

    return new_path
#---------------------------------------
def _test_add_new_dir() -> bool:
    dirs = ("testdir1", "testdir2", "testdir3")
    home = os.path.expanduser("~")
    try:
        new_path = add_new_dir(home, *dirs)
        return os.path.isdir(new_path)
    finally:
        os.system("rm -r ~/testdir1")


########################################


class BadYamlWarning(UserWarning):
    """Class of warnings related to .yaml files reading."""


########################################


def yaml_loader(file_path: PathLike, on_error: Literal["raise", "warn", "ignore"] = "warn"):
    """Load a YAML file and casts it explicitly to a dict"""
    if not isinstance(file_path, (Path, str)):
        raise TypeError(f"Expected a Path object or str, got {type(file_path)} instead.")

    assert str(file_path).endswith(".yaml"), \
    f"{str(file_path)} is not a .yaml file format."

    file_path = Path(file_path)

    assert file_path.is_file(), \
    f"{str(file_path)}: no such file found."

    yaml = YAML()
    with open(file_path, encoding="utf-8") as yaml_file:
        try:
            yaml_data = yaml.load(yaml_file)
        except Exception as exc:
            if on_error == "raise":
                raise exc
            if on_error == "warn":
                warnings.warn(
                    f"An exception was thrown during yaml loading of file {str(file_path)}.\n"
                    f"Data written in this file is ignored to proceed.\n"
                    f"Thrown exception below:\n"
                    f"{exc}", BadYamlWarning
                )
            return {}
        return dict(yaml_data)


########################################

if __name__ == "__main__":
    if _test_add_new_dir():
        print("add_new_dir test passed !")