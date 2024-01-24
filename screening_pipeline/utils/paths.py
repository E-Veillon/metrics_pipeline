'''
Implements functions to manage and operate on paths.
'''

from typing import Union, Sequence, List
from pathlib import Path
from itertools import repeat

PathLike = Union[str, Path]

def add_new_dir(base_dir: PathLike, new_dir_name: PathLike) -> Path:
    '''
    Creates a new directory inside given base directory.
    
    Parameters:
        base_dir (str|Path): The base directory inside which the new one will be created.

        new_dire_name (str): The name of the new subdirectory to create.

    Returns:
        Path: The Path object pointing to the new subdirectory.
    '''

    assert isinstance(base_dir, PathLike)
    assert Path(base_dir).is_dir()
    assert isinstance(new_dir_name, PathLike)

    new_dir  = Path('/'.join(str(base_dir), new_dir_name))
    new_dir.mkdir()

    return new_dir

def batch_add_new_dirs(
        base_dir: PathLike, 
        new_subdirs: Sequence[PathLike]
    ) -> List[Path]:
    '''
    Iterates through new_subdirs to create a bunch of subdirectories in base_dir.
    Returns the list of created paths.

    Parameters:
        base_dir (str|Path):    An existing directory where subdirectories should be created.

        new_subdirs (str|Path): Names or subpaths relative to base_dir to create directories at.
    
    Returns:
        List[Path]: A list of all newly created paths starting from base_dir.
    '''

    assert isinstance(base_dir, PathLike)

    base_dir = Path(base_dir)
    
    assert base_dir.is_dir()
    assert all(isinstance(subdir, PathLike) for subdir in new_subdirs)

    new_dirs = list(map(add_new_dir, repeat(base_dir), new_subdirs))
    
    return new_dirs