"""I/O operations with YAML files."""

import typing as tp
import warnings

from ruamel.yaml import YAML

from .io_base import PathLike, check_file_or_dir


class BadYamlWarning(UserWarning):
    """Class of warnings related to .yaml files reading."""


def load_yaml_as_dict(
    file_path: PathLike, on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
) -> dict:
    """Load a YAML file and casts it explicitly to a dict."""
    check_file_or_dir(file_path, "file",  allowed_formats="yaml")

    yaml = YAML()
    with open(file_path, "rt", encoding="utf-8") as yaml_file:
        try:
            yaml_data = yaml.load(yaml_file)
        except Exception as exc:
            if on_error == "raise":
                raise exc
            if on_error == "warn":
                warnings.warn(
                    f"An exception was thrown during yaml loading of file {str(file_path)}. "
                    "Data written in this file is ignored to proceed.\nThrown exception below:\n"
                    f"{exc}",
                    BadYamlWarning
                )
            return {}
        return dict(yaml_data)