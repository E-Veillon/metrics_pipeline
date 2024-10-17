#!/usr/bin/python
"""
Module implementing general utilitary functions.
"""


import itertools
from typing import Sequence, List, Any


########################################


def flatten(sequence: Sequence[Any], level_of_flattening: int = 1) -> List:
    """
    Unpacks a nested sequence without modifying elements order.

    Parameters:
        sequence (Sequence):        The iterable to unpack.

        level_of_flattening (Int):  The number of nested levels to unpack.
                                    Defaults to 1.
    
    Returns:
        A list flattened the specified number of times.
    """
    for _ in range(1, level_of_flattening + 1):
        sequence = list(itertools.chain.from_iterable(sequence))
    return sequence
