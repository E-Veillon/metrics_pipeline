"""
A minimalist way to visualize progression while iterating through long sequences.
Possibility to custom printed messages.
"""


import typing as tp
import typing_extensions as tpe
import warnings


class VisualIterator:
    """
    A minimalist way to visualize progression while iterating through long sequences.
    Possibility to custom printed messages.
    """
    def __init__(
        self: tpe.Self,
        iterable: tp.Iterable[tp.Any],
        desc: str | None = None,
        unit: str | None = None,
        end_desc: str | None = None,
        percent: bool = False,
        significant_figures: int = 2,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> None:
        """
        A minimalist way to visualize progression while iterating through long sequences.

        Parameters
        ----------

        iterable: Iterable
            An iterable object. If you need to wrap a lazy iterator or generator,
            do not use the constructor directly but appropriate classmethod instead.

        desc: str
            Description to put before the counter that describes the task.
            Defaults to "Iterating".

        unit: str, optional
            name of the iterated objects. Defaults to "iterations".

        end_desc: str, optional
            Message to print when the task is finished. Defaults to "Done.".

        percent: bool
            Whether to print a percentage of progression next to the counter.
            Defaults to False.

        significant_figures: int
            Number of decimal places to show in the percentage (if activated).
            Defaults to 2.

        on_error: str
            What to do in case of an exception raising while trying to wrap the iterable.
            Defaults to "warn".
        """
        self.iterable = iterable
        self.desc = "Iterating" if desc is None else desc
        self.unit = "iterations" if unit is None else unit
        self.end_desc = "Done." if end_desc is None else end_desc
        self.percent = percent
        self.significant_figures = significant_figures
        self.on_error: tp.Literal["raise", "warn", "ignore"] = on_error
        self.default_length = -1
        self._length = self.default_length

        try:
            self.length = len(iterable) # type: ignore
        except TypeError:
            err_msg = (
                "Could not access the length of given iterable object. "
                "If you are trying to wrap a lazy iterator or generator object, "
                "make sure to use the appropriate classmethod and not the "
                f"constructor directly."
            )
            warn_msg = err_msg + (
                f" Length is now set back to its default value of {self.default_length}, "
                "and percentage printing will be deactivated."
            )
            raise_or_warn(self.on_error, TypeError, err_msg, UserWarning, warn_msg)
            self.length = self.default_length

    @property
    def length(self) -> int:
        return self._length
    
    @length.setter
    def length(self, value: int) -> None:
        match value:
            case int() if value >= -1:
                self._length = value

            case int():
                err_msg = (
                    f"{self.__class__.__name__}: "
                    f"Length cannot be inferior to default of {self.default_length}, "
                    "used for undetermined length."
                )
                warn_msg = err_msg + (
                    f" Setting length back to default value of {self.default_length}."
                )
                raise_or_warn(self.on_error, ValueError, err_msg, UserWarning, warn_msg)
                self._length = -1

            case _:
                err_msg = (
                    f"{self.__class__.__name__}: "
                    "'length' expected a type 'int', "
                    f"got '{type(value).__name__}' instead."
                )
                warn_msg = err_msg + (
                    f" Setting length back to default value of {self.default_length}."
                )
                raise_or_warn(self.on_error, TypeError, err_msg, UserWarning, warn_msg)
                self._length = -1

    @length.deleter
    def length(self) -> None:
        self._length = self.default_length

    def __getattr__(self, attr: str) -> tp.Any:
        try:
            val = getattr(self.iterable, attr)
        except AttributeError as exc:
            raise AttributeError(
                f"Attribute {attr!r} is not defined, neither for the "
                f"{self.__class__.__name__!r} object nor for the iterable inside it."
            ) from exc
        else:
            return val

    def __len__(self: tpe.Self) -> int:
        return self.length

    def __iter__(self: tpe.Self) -> tp.Any:
        for idx, obj in enumerate(self.iterable, start=1):
            maxlen = "Unknown" if self.length < 0 else self.length
            counter = f"{self.desc}: {idx}/{maxlen} {self.unit}"

            if self.percent and isinstance(maxlen, int):
                sf = self.significant_figures
                fmt = f".{sf}f"
                percent = f" ({round(idx / maxlen * 100, sf):{fmt}}%)"
                print(counter + percent, end="\r", flush=True)

            else:
                print(counter, end="\r", flush=True)

            yield obj

        print(f"\n{self.end_desc}")

    def __repr__(self: tpe.Self) -> str:
        return str(self) + f"\nWrapped iterable: {self.iterable}"
    
    def __str__(self: tpe.Self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"length={self.length},desc={self.desc},"
            f"end_desc={self.end_desc},percent={self.percent},"
            f"significant_figures={self.significant_figures},"
            f"on_error={self.on_error})"
        )

    @classmethod
    def from_iterator(
        cls, iterator: tp.Iterator|tp.Generator, **kwargs
    ) -> tpe.Self:
        """
        Get a VisualIterator wrapping a lazy iterator instead of an explicit sequence.

        WARNING: This method converts the iterator into list in order to get its actual
        length. If the iterator to wrap can generate enough data to overflow memory,
        use from_big_iterator instead for lazy wrapping.

        Parameters:
            iterator (Iterator):    The iterator to wrap.

            **kwargs:               Keyword arguments to pass to the constructor.

        Returns:
            VisualIterator.
        """
        return cls(list(iterator), **kwargs)

    @classmethod
    def from_big_iterator(
        cls, iterator: tp.Iterator|tp.Generator, n_elts: int | None = None, **kwargs
    ) -> tpe.Self:
        """
        Get a VisualIterator wrapping a lazy iterator instead of an explicit sequence.

        This method wraps iterators lazily, so that they do not overflow memory.
        
        Parameters:
            iterator (Iterator):    The iterator to wrap.
            
            n_elts (int):           The total number of elements to consider to render
                                    a max limit and a percentage. This method won't search
                                    for a total number of elements in the iterator as it
                                    would overflow memory, hence it will print a maximum of
                                    '-1' and deactivate the 'percent' option if no value
                                    is given.
            
            **kwargs:               Keyword arguments to pass to the constructor.

        Returns:
            VisualIterator.    
        """
        kwargs.setdefault("on_error", "ignore")
        visual_iterator = cls(iterator, **kwargs)

        if isinstance(n_elts, int) and n_elts > 0:
            visual_iterator.length = n_elts

        return visual_iterator


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

    exception: Exception
        What class of exception to raise when `raise_or_warn` is set to "raise".
        If not given and `raise_or_warn` is set to "raise", defaults to base
        exception class `Exception`.

    err_msg: str | None
        The message to print when raising an exception.

    warn_type: Warning
        What class of warning to show when `raise_or_warn` is set to "warn".
        If not given and `raise_or_warn` is set to "xarn", defaults to
        warning class `UserWarning`.

    warn_msg: str | None
        The message to print when printing a warning.
    
    NOTE: if the message argument corresponding to given `raise_or_warn` action is not given,
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
        warn_msg = warn_msg if warn_msg is not None else err_msg
        warnings.warn(warn_msg, warn_type, stacklevel = 2)


def test_visual_iterator() -> None:
    visual_construct = VisualIterator(list(range(5)), percent=True)
    visual_from_iter = VisualIterator.from_iterator(iter(range(5)), percent=True)
    visual_from_big_iter = VisualIterator.from_big_iterator(iter(range(100)), n_elts=100, percent=True)

    for vit in (visual_construct, visual_from_iter, visual_from_big_iter):
        if vit is visual_construct or vit is visual_from_iter:
            assert vit.iterable == [0, 1, 2, 3, 4]
            assert vit.length == 5

        if vit is visual_from_big_iter:
            assert vit.iterable == iter(range(100))
            assert vit.length == 100

        assert len(vit) == vit.length
        assert vit.default_length == -1
        assert vit.desc == "Iterating"
        assert vit.end_desc == "Done."
        assert vit.percent == True
        assert vit.significant_figures == 2
        assert vit.on_error == "warn"

