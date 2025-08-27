#!/usr/bin/python
"""This module defines common guard clauses used throughout the pipeline."""


import typing as typ


CompareStr = typ.Union[
    typ.Literal["=="],
    typ.Literal["!="],
    typ.Literal["<="],
    typ.Literal[">="],
    typ.Literal["<"],
    typ.Literal[">"],
]


def _get_class_name(class_str: str) -> str:
    """Extract the class name from a string of the form '<class 'name'>'."""
    return class_str.split()[-1].rstrip(">")


def check_type(
    obj: typ.Any, obj_name: str, wanted_types: typ.Tuple[typ.Type, ...]
) -> None:
    """Check whether the type of passed object matches wanted type."""
    if not isinstance(obj, wanted_types):
        obj_type = _get_class_name(str(type(obj)))

        wanted_str_types = ""
        for allowed_type in wanted_types[:-1]:
            str_type = _get_class_name(str(allowed_type))
            wanted_str_types += f"{str_type}, "
        wanted_str_types += f"or {_get_class_name(str(wanted_types[-1]))}"

        raise TypeError(
            f"'{obj_name}' argument expected a type {wanted_str_types}, "
            f"got {obj_type} instead."
        )

def check_num_value(
    val: int|float|None, val_name: str, cdt: CompareStr = "==", ref_val: int|float = 0
) -> None:
    """
    Test given condition on given numeric value (int or float),
    and raises a preformatted ValueError when the condition is not met.

    Parameters:
        val (int|float):        Numerical value to check.

        val_name (str):         name of the variable passed to 'val'.

        cdt (str):              String representing a numerical comparison
                                operator (e.g. ">="). Defaults to "==".

        ref_val (int|float):    The value to compare 'val' to. Defaults to 0.
    """
    check_type(val, "val", (int, float))
    assert val is not None, "Never triggered, used for type checker."
    check_type(val_name, "val_name", (str,))
    check_type(ref_val, "ref_val", (int, float))

    match cdt:
        case "==":
            if val != ref_val:
                raise ValueError(
                    f"'{val_name}' should be equal to {ref_val}, "
                    f"got {val}."
                )
        case "!=":
            if val == ref_val:
                raise ValueError(
                    f"'{val_name}' should be different from {ref_val}, "
                    f"got {val}."
                )
        case "<=":
            if val > ref_val:
                raise ValueError(
                    f"'{val_name}' should be inferior or equal to {ref_val}, "
                    f"got {val}."
                )
        case ">=":
            if val < ref_val:
                raise ValueError(
                    f"'{val_name}' should be superior or equal to {ref_val}, "
                    f"got {val}."
                )
        case ">":
            if val <= ref_val:
                raise ValueError(
                    f"'{val_name}' should be strictly superior to {ref_val}, "
                    f"got {val}."
                )
        case "<":
            if val >= ref_val:
                raise ValueError(
                    f"'{val_name}' should be strictly inferior to {ref_val}, "
                    f"got {val}."
                )
        case str():
            raise NotImplementedError(
                f"'{cdt}' is not a supported comparison operator."
            )
        case _:
            check_type(cdt, "cdt", (str,))


if __name__ == "__main__":
    NUM = 5
    FLOAT = 3.14
    BOOL = True
    STR = "hello"
    SET = {1, 2, 3}
    TUPLE = ("cdo", 2)
    LIST = [1, 3.2, "hello"]
    DICT = {"a": 1, "b": 2, "c": 3}
    check_type(NUM, "num", (int,))
    check_type(FLOAT, "float_val", (float, int))
    check_type(BOOL, "bol", (bool,))
    check_type(STR, "string", (str,))
    check_type(SET, "sett", (typ.Set,))
    check_type(SET, "sett", (set,))
    check_type(TUPLE ,"tup", (typ.Tuple,))
    check_type(TUPLE, "tup", (tuple,))
    check_type(LIST, "lst", (typ.List,))
    check_type(LIST, "lst", (list,))
    check_type(DICT, "dct", (typ.Dict,))
    check_type(DICT, "dct", (dict,))
    check_num_value(NUM, "num", "==", 5)
    check_num_value(FLOAT, "float_val", "<", 4.2)
    try:
        check_num_value(NUM, "num", ">", 5)
    except ValueError as exc:
        print("Wrong check_num_value test successfully failed for 'num'!")
    try:
        check_num_value(FLOAT, "float_val", ">=", 3.5)
    except ValueError as exc:
        print("Wrong check_num_value test successfully failed for 'float_val'!")
    print("All tests passed !")
