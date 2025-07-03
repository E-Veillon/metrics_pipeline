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
        end_desc: str | None = None,
        percent: bool = False,
        significant_figures: int = 2,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> None:
        """
        Parameters:
            iterable (Iterable):        An iterable object. If you need to wrap a lazy iterator
                                        or generator, do not use the constructor directly but
                                        appropriate classmethod instead.

            desc (str):                 Description to put before the counter that describes the
                                        task. Defaults to "Iterating".

            end_desc (str):             Message to print when the task is finished.
                                        Defaults to "Done.".

            percent (bool):             Whether to print a percentage of progression next to the
                                        counter. Defaults to False.

            significant_figures (int):  Number of decimal places to show in the percentage
                                        (if activated). Defaults to 2.

            on_error (str):             What to do in case of an exception raising while trying
                                        to wrap the iterable. Defaults to "warn".
        """
        self.iterable = iterable
        self.desc = desc if desc is not None else "Iterating"
        self.end_desc = end_desc if end_desc is not None else "Done."
        self.percent = percent
        self.significant_figures = significant_figures
        self.on_error = on_error
        self._length = -1 # Default value

        try:
            self.length = len(iterable) # type: ignore
        except TypeError as exc:
            self.length = -1
            self.percent = False
            error_msg = (
                "Could not access the length of given iterable object. "
                "If you are trying to wrap a lazy iterator or generator object, "
                "make sure to use the appropriate classmethod and not the "
                f"constructor directly. Default length is set to {self.length}, "
                "and percentage printing is deactivated as it will make no sense."
            )
            if self.on_error == "raise":
                raise TypeError(error_msg) from exc
            if self.on_error == "warn":
                warnings.warn(error_msg)

    @property
    def length(self) -> int:
        return self._length
    
    @length.setter
    def length(self, value: int) -> None:
        if isinstance(value, int) and value >= -1:
            self._length = value

    @length.deleter
    def length(self) -> None:
        self._length = -1

    def __len__(self: tpe.Self) -> int:
        return self.length

    def __iter__(self: tpe.Self) -> tp.Any:
        for idx, obj in enumerate(self.iterable, start=1):
            counter = f"{self.desc}: {idx}/{self.length}"

            if self.percent:
                sf = self.significant_figures
                fmt = f".{sf}f"
                percent = f" ({round(idx / self.length * 100, sf):{fmt}}%)"
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
