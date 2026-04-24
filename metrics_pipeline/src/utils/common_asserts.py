#!/usr/bin/python
"""This module defines common guard clauses used throughout the pipeline."""


import typing as tp
import warnings

CompareStr = tp.Union[
    tp.Literal["=="],
    tp.Literal["!="],
    tp.Literal["<="],
    tp.Literal[">="],
    tp.Literal["<"],
    tp.Literal[">"],
]

def check_type(
    obj: tp.Any, obj_name: str, wanted_types: tp.Tuple[tp.Type, ...]
) -> None:
    """Check whether the type of passed object matches wanted type."""
    if isinstance(obj, wanted_types):
        return

    obj_type = type(obj).__name__

    wanted_str_types = ""
    for allowed_type in wanted_types[:-1]:
        str_type = allowed_type.__name__
        wanted_str_types += f"{str_type}, "
    wanted_str_types += f"or {wanted_types[-1].__name__}"

    raise TypeError(
        f"'{obj_name}' argument expected a type {wanted_str_types}, "
        f"got {obj_type} instead."
    )

def check_num_value(
    val: int | float, val_name: str, cdt: CompareStr = "==", ref_val: int | float = 0
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
    check_type(val, val_name, (int, float))
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


def raise_or_warn(
    action: tp.Literal["raise", "warn", "ignore"],
    exception: type[Exception] | None = None,
    err_msg: str | None = None,
    warn_type: type[Warning] | None = None,
    warn_msg: str | None = None
) -> None:
    """
    Flexibly raise an exception, print a warning or ignore
    with a message that can be different in each case.

    Parameters
    ----------

    action: "raise" | "warn" | "ignore"
        What kind of action to do when going through this function.

    exception: Exception, optional
        What class of exception to raise when `raise_or_warn` is set to "raise".
        If not given and `raise_or_warn` is set to "raise", defaults to base
        exception class `Exception`.

    err_msg: str, optional
        The message to print when raising an exception.

    warn_type: Warning, optional
        What class of warning to show when `raise_or_warn` is set to "warn".
        If not given and `raise_or_warn` is set to "warn", defaults to
        warning class `UserWarning`.

    warn_msg: str, optional
        The message to print when printing a warning.

    Notes
    -----
    If the message argument corresponding to given `raise_or_warn` action is not given,
    the message set for the other action is used instead, so if the message is the same for
    both actions it can be passed only once to one or the other message argument indifferently.
    At least one message argument must be given.
    """
    match (err_msg, warn_msg):
        case (str(), str()): pass
        case (None, str()): err_msg = warn_msg
        case (str(), None): warn_msg = err_msg
        case (None, None):
            raise ValueError(
                f"{raise_or_warn.__name__}: at least one of either 'err_msg' or 'warn_msg' "
                "arguments must be set."
            )
        case _:
            raise TypeError(
                f"{raise_or_warn.__name__}: At least one of either 'err_msg' or 'warn_msg' "
                "arguments were given a wrong type:\n"
                f"- 'err_msg' expected a type 'str', got {type(err_msg).__name__!r}.\n"
                f"- 'warn_msg' expected a type 'str', got {type(warn_msg).__name__!r}.\n"
            )

    if action == "raise":
        exception = exception if exception is not None else Exception
        err_msg = err_msg if err_msg is not None else warn_msg
        raise exception(err_msg)

    if action == "warn":
        warn_type = warn_type if warn_type is not None else UserWarning
        warn_msg = warn_msg if warn_msg is not None else err_msg
        warnings.warn(warn_msg, warn_type, stacklevel = 2)


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
    check_type(SET, "sett", (tp.Set,))
    check_type(SET, "sett", (set,))
    check_type(TUPLE ,"tup", (tp.Tuple,))
    check_type(TUPLE, "tup", (tuple,))
    check_type(LIST, "lst", (tp.List,))
    check_type(LIST, "lst", (list,))
    check_type(DICT, "dct", (tp.Dict,))
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
