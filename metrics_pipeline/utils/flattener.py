#!/usr/bin/python
"""
Module implementing general utilitary functions.
"""


import itertools as it
import typing as tp


########################################


def flatten(
    sequence: tp.Sequence[tp.Any], level_of_flattening: int = 1
) -> list[tp.Any]:
    """
    Unpacks a nested sequence without modifying elements order.

    Parameters:
        sequence (Sequence):        The iterable to unpack.

        level_of_flattening (Int):  The number of nested levels to unpack.
                                    Defaults to 1.
    
    Returns:
        A list flattened the specified number of times.
    """
    result = list(sequence)
    for _ in range(1, level_of_flattening + 1):
        result = list(it.chain.from_iterable(result))
    return result
