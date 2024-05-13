'''
Module implementing general utilitary functions.
'''


from pathlib import Path
from typing import Sequence, List, Any, Literal
import itertools
from ruamel.yaml import YAML
import warnings

from screening_pipeline.utils.custom_types import PathLike

class BadYamlWarning(UserWarning):
    '''Class of warnings related to .yaml files reading.'''

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


def _yaml_loader(file_path: PathLike, on_error: Literal['raise', 'warn', 'ignore'] = 'warn'):
    if not isinstance(file_path, (Path, str)):
        raise TypeError(f"Expected a Path object or str, got {type(file_path)} instead.")
    
    assert str(file_path).endswith('.yaml'), \
    f'{str(file_path)} is not a .yaml file format.'

    file_path = Path(file_path)

    assert file_path.is_file(), \
    f'{str(file_path)}: no such file found.'

    yaml = YAML()
    with open(file_path, encoding="utf-8") as yaml_file:
        try: yaml_data = yaml.load(yaml_file)
        except Exception as exc:
            if on_error == 'raise': raise exc
            elif on_error == 'warn':
                warnings.warn(
                    f'An exception was thrown during yaml loading of file {str(file_path)}.\n'
                    f'Data written in this file is ignored to proceed.\n'
                    f'Thrown exception below:\n'
                    f'{exc}', BadYamlWarning
                )
            return {}
        return dict(yaml_data)

#def _cast_str_to_seq(string: str, /, *, seq: Literal['tuple','list']) -> Union[Tuple, List]:
#   string = string[1:-1]
#   if string.find('[') != -1 or string.find('(') != -1 or string.find('{') != -1:
#       #TODO: gérer les conteneurs imbriquées
#       raise NotImplementedError
#   else:
#       new_seq = tuple(string.split(sep=',')) if seq == 'tuple' else string.split(sep=',')
#       return new_seq
#
#def _parse_dict(dct: Dict[str, str]) -> Dict[str, Any]:
#   
#   for k, v in dct.items():
#       if isinstance(v, str):
#           if v.lower() == 'true': dct[k] = True
#           elif v.lower() == 'false': dct[k] = False
#           elif v.isdecimal(): dct[k] = int(v)
#           elif is_float(v): dct[k] = float(v)
#           elif v.startswith('[') and v.endswith(']'):
#               v = _cast_str_to_seq(v, seq='list')
#               dct[k] = _parse_seq(v, seq='list')
#           elif v.startswith('(') and v.endswith(')'):
#               v = _cast_str_to_seq(v, seq='tuple')
#               dct[k] = _parse_seq(v, seq='tuple')
#           elif v.startswith('{') and v.endswith('}') and ':' in v:
#               v = _cast_str_to_dict(v)
#               dct[k] = _parse_dict(v)
#           elif v.startswith('{') and v.endswith('}'):
#               v = _cast_str_to_set(v)
#               dct[k] = _parse_set(v)




