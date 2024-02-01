'''
Module implementing general utilitary functions.
'''


from pathlib import Path
from typing import Sequence, List, Any
import itertools
from ruamel.yaml import YAML

from screening_pipeline.utils.typing import PathLike

def is_float(string: str) -> bool:
    str_list = string.split(sep='.')
    return ('.' in string) and (len(str_list) <= 2) and all([nbr.isdecimal() for nbr in str_list])


def flatten(sequence: Sequence[Any], level_of_flattening: int = 1) -> List:
    '''
    Unpacks a nested sequence without modifying elements order.

    Parameters:
        sequence (Sequence):        The iterable to unpack.

        level_of_flattening (Int):  The number of nested levels to unpack.
                                    Defaults to 1.
    
    Returns:
        A list flattened the specified number of times.
    '''

    for _ in range(1, level_of_flattening + 1):
        sequence = list(itertools.chain.from_iterable(sequence))
    return sequence


def _yaml_loader(file_path: PathLike):
    if not isinstance(file_path, (Path, str)):
        raise TypeError(f"Expected a Path object or str, got {type(file_path)} instead.")
    
    assert str(file_path).endswith('.yaml'), \
    f'{str(file_path)} is not a .yaml file format.'

    file_path = Path(file_path)

    assert file_path.is_file(), \
    f'{str(file_path)}: no such file found.'

    yaml = YAML()
    with open(file_path, encoding="utf-8") as yaml_file:
        try: yaml_data = yaml.load(yaml_file) or {}
        except Exception as exc:
            warn_msg = f'An exception was thrown during yaml loading of file {str(file_path)}.\n\
                        Data written in this file is ignored to proceed.\n\
                        Thrown exception below:\n\
                        {exc}'
            print(warn_msg)
            return {}
        return dict(yaml_data)

