'''
Implements functions to manage and operate on paths.
'''

from typing import Union, Sequence, List
from pathlib import Path
from functools import partial

PathLike = Union[str, Path]

def add_new_dir(base_dir: PathLike, new_dir_name: str) -> Path:
    '''
    Creates a new directory inside given base directory.
    
    Parameters:
        base_dir (str|Path): The base directory inside which the new one will be created.
        new_dire_name (str): The name of the new subdirectory to create.

    Returns:
        Path: The Path object pointing to the new subdirectory.
    '''

    assert isinstance(base_dir, PathLike)

    base_dir = Path(base_dir)
    
    assert base_dir.exists() and base_dir.is_dir()
    assert isinstance(new_dir_name, str)

    new_dir  = base_dir / new_dir_name
    new_dir.mkdir()

    return new_dir

def batch_add_new_dirs(
        base_dir: PathLike, 
        new_subdirs: Sequence[str]
    ) -> List[PathLike]:
    '''
    Iterates through new_subdirs to create a bunch of subdirectories in base_dir.
    Returns the list of created paths.
    '''

    assert isinstance(base_dir, PathLike)

    base_dir = Path(base_dir)
    
    assert base_dir.exists() and base_dir.is_dir()
    assert all(isinstance(subdir, str) for subdir in new_subdirs)
    
    based_add_new_dir = partial(add_new_dir, base_dir=base_dir)
    new_dirs = list(map(based_add_new_dir, new_subdirs))
    
    return new_dirs