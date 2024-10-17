#!/usr/bin/python
"""
Implements functions to manage and operate on paths.
"""


import os
import warnings
from typing import Literal, Sequence, List, Tuple
from pathlib import Path
from ruamel.yaml import YAML

# LOCAL IMPORTS
from .common_asserts import check_type
from .custom_types import PathLike

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


def check_file_format(filename: PathLike, *, allowed_formats: str|Sequence[str]) -> None:
    """
    Verify that extension format of given file and wanted file format match,
    no matter if given file exists or not.

    Parameters:
        filename (str|Path):            Name or path of the file to check.

        allowed_formats (str|[str]):    Allowed extension formats for the checked file,
                                        without the dot separator (e.g. "txt" and not ".txt").
                                        Can be given in a single string or as a sequence
                                        (list or tuple) of strings.
    """
    check_type(filename, "filename", (str, Path))
    check_type(allowed_formats, "allowed_formats", (str, Sequence))
    if isinstance(allowed_formats, str):
        allowed_formats = (allowed_formats,)
    else:
        for idx, ext in enumerate(allowed_formats):
            check_type(ext, f"structures[{idx}]", (str,))

    filename = str(filename)
    file_format = filename.rsplit(sep=".", maxsplit=1)[-1]

    if all(file_format != ext for ext in allowed_formats):
        plural = "s are" if len(allowed_formats) > 1 else " is"
        formats_str = ", ".join(["'" + ext + "'" for ext in allowed_formats])
        raise ValueError(
            f"{filename}: allowed file format{plural} {formats_str}, "
            f"got '{file_format}' format instead."
        )


########################################


def check_file_or_dir(
    path: PathLike,
    file_or_dir: Literal["file", "dir"] = "file",
    *,
    allowed_formats: str|Sequence|None = None
) -> None:
    """
    Verify existence and optionally extension format of given path.
    
    Parameters:
        path (str|Path):                Path to verify.

        file_or_dir (str):              Whether the path should lead to a file or a directory.
                                        If the path exists but is not the right data type,
                                        an error will still be raised for not finding it.

        allowed_formats (str|[str]):    If the path should lead to a file with a specific format
                                        extension, provide here wanted extension without the dot
                                        separator (e.g. "txt" and not ".txt"). If several formats
                                        are possible, give a sequence (list or tuple) of them.
    """
    check_type(path, "path", (str, Path))
    assert file_or_dir in {"file", "dir"}
    if allowed_formats is not None:
        check_type(allowed_formats, "allowed_formats", (str, List, Tuple))

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
        and allowed_formats is not None
    ):
        check_file_format(path, allowed_formats=allowed_formats)


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
    check_file_or_dir(base_dir, "dir")
    for idx, new_dir in enumerate(new_dirs):
        check_type(new_dir, f"new_dirss[{idx}]", (str,))

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
    check_file_or_dir(file_path, "file",  allowed_formats="yaml")

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
