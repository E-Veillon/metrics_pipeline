"""Basic I/O operations and validations."""

import os
from pathlib import Path
import typing as tp


PathLike = Path | str

# Main paths inside pipeline file tree
ROOT = Path(__file__).resolve().parent.parent.parent
"""Absolute path to the GenMat main directory."""

SRCPATH = ROOT / "src"
"""Absolute path to the GenMat src directory."""

CONFIGPATH = ROOT / "config"
"""Absolute path to the GenMat config directory."""


class EmptyDirectoryError(FileNotFoundError):
    """Directory is empty."""

def check_file_format(
    filename: PathLike | None, *, allowed_formats: str | tuple[str, ...]
) -> None:
    """
    Verify that extension format of given file and wanted file format match,
    no matter if given file exists or not.

    Parameters
    ----------
    filename: str | Path
        Name or path of the file to check.

    allowed_formats: str | tuple[str]
        Allowed extension formats for the checked file, without the dot separator
        (e.g. "txt" and not ".txt").
    """
    assert isinstance(filename, (str, Path)), TypeError(
        f"'filename' expected a type 'str' or 'Path', got {type(filename).__name__!r}."
    )
    assert all(isinstance(fmt, str) for fmt in allowed_formats), TypeError(
        "'allowed_formats' expected a type 'str' or 'tuple[str]', "
        f"got {type(filename).__name__!r}."
    )
    if isinstance(allowed_formats, str):
        allowed_formats = (allowed_formats,)

    filename = str(filename)
    file_format = filename.rsplit(sep=".", maxsplit=1)[-1]

    if all(file_format != ext for ext in allowed_formats):
        plural = "s are" if len(allowed_formats) > 1 else " is"
        formats_str = ", ".join(["'" + ext + "'" for ext in allowed_formats])
        raise ValueError(
            f"{filename}: allowed file format{plural} {formats_str}, "
            f"got '{file_format}' format instead."
        )


def check_file_or_dir(
    path: PathLike | None,
    file_or_dir: tp.Literal["file", "dir"] = "file",
    *,
    check_empty: bool = False,
    allowed_formats: str| tuple[str, ...] | None = None
) -> None:
    """
    Verify existence and optionally extension format of given path.
    
    Parameters
    ----------
        path: str | Path
            Path to verify.

        file_or_dir: str
            Whether the path should lead to a file or a directory.
            If the path exists but is not the right data type, `FileNotFoundError` is raised.

        check_empty: bool
            Whether to check if the directory is empty. Raises `EmptyDirectoryError` if it is
            the case. Ignored for file checking. Defaults to False.

        allowed_formats: str | tuple[str]
            If the path should lead to a file with a specific format extension, provide here
            wanted extension without the dot separator (e.g. "txt" and not ".txt").
            If several formats are possible, give a tuple of them.
    """
    assert isinstance(path, (str, Path)), TypeError(
        f"'filename' expected a type 'str' or 'Path', got {type(path).__name__!r}."
    )
    assert file_or_dir in {"file", "dir"}, ValueError(
        f"'file_or_dir' only supports 'file' or 'dir', got {file_or_dir!r}."
    )
    path = str(path)

    if file_or_dir == "dir" and not os.path.isdir(path):
        raise FileNotFoundError(
            f"{path}: No such directory found."
        )
    if file_or_dir == "dir" and check_empty and not os.listdir(path):
        raise EmptyDirectoryError(
            f"{path}: Directory exists but is empty."
        )
    if file_or_dir == "file" and not os.path.isfile(path):
        raise FileNotFoundError(
            f"{path}: No such file found."
        )
    if file_or_dir == "file" and allowed_formats is not None:
        check_file_format(path, allowed_formats=allowed_formats)