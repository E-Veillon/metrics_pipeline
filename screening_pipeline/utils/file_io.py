"""
Implements functions to manage and operate on paths.
"""


import os
import warnings
from itertools import repeat
from typing import Sequence, List, Literal
from pathlib import Path
from ruamel.yaml import YAML


from screening_pipeline.utils import PathLike


########################################


def add_new_dir(base_dir: PathLike, new_dir_name: PathLike) -> Path:
    """
    Creates a new directory inside given base directory. If the subdirectory already exists, 
    it is untouched but its path is still returned.
    
    Parameters:
        base_dir (str|Path): The base directory inside which the new one will be created.

        new_dire_name (str): The name of the new subdirectory to create.

    Returns:
        Path: The Path object pointing to the new subdirectory.
    """

    assert isinstance(base_dir, PathLike)
    assert os.path.isdir(str(base_dir))
    assert isinstance(new_dir_name, PathLike)

    new_dir = os.path.join(str(base_dir), str(new_dir_name))
    os.makedirs(new_dir, exist_ok=True)

    return new_dir

########################################

def batch_add_new_dirs(
        base_dir: PathLike, 
        new_subdirs: Sequence[PathLike]
    ) -> List[Path]:
    """
    Iterates through new_subdirs to create a bunch of subdirectories in base_dir.
    Returns the list of created paths. If the subdirectory already exists, it is 
    untouched but its path is still returned.

    Parameters:
        base_dir (str|Path):    An existing directory where subdirectories should be created.

        new_subdirs (str|Path): Names or subpaths relative to base_dir to create directories at.
    
    Returns:
        List[Path]: A list of all newly created paths starting from base_dir.
    """

    assert isinstance(base_dir, PathLike)

    base_dir = Path(base_dir)
    
    assert base_dir.is_dir()
    assert all(isinstance(subdir, PathLike) for subdir in new_subdirs)

    new_dirs = list(map(add_new_dir, repeat(base_dir), new_subdirs))
    
    return new_dirs


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
            elif on_error == "warn":
                warnings.warn(
                    f"An exception was thrown during yaml loading of file {str(file_path)}.\n"
                    f"Data written in this file is ignored to proceed.\n"
                    f"Thrown exception below:\n"
                    f"{exc}", BadYamlWarning
                )
            return {}
        return dict(yaml_data)


########################################