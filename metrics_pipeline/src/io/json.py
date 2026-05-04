"""Load and write JSON files with extra features."""

import typing as tp
import json
from pathlib import Path

import numpy as np

from .io_base import check_file_or_dir, check_file_format, PathLike


class JsonLoader:
    """Load JSON files with extra features."""
    def __init__(self, filepath: PathLike) -> None:
        """
        Load JSON files with extra features.
        
        Parameters
        ----------
        filepath: str | Path
            Path to the file to load.
        """
        check_file_or_dir(filepath, "file", allowed_formats="json")
        self.filepath = Path(filepath)

    def load(self) -> tp.Any:
        """Load the JSON data as-is."""
        with self.filepath.open("rt", encoding="utf-8") as fp:
            data = json.load(fp)
        return data

    def load_as_dict(self) -> dict[tp.Any, tp.Any]:
        """
        Load the file and cast the data into a dict explicitly.
        - If it is already a dict, it will be returned untouched.
        - If it is a list, it will be cast into a dict with integer
        position indices as integer keys.
        - In any other case, the whole data is stored in a one-element dict
        at key 0 (integer).
        """
        data = self.load()

        if isinstance(data, dict):
            return data

        elif isinstance(data, list):
            return {k: v for k, v in enumerate(data)}

        else:
            return {0: data}

    def load_as_list(self) -> list[tp.Any]:
        """
        Load the file and cast the data into a list explicitly.
        - If it is already a list, it is returned untouched.
        - If it is a dict, it is cast into a list of key-value pairs tuples.
        - In any other case, all the data is stored in a one-element list.
        """
        data = self.load()

        if isinstance(data, list):
            return data

        if isinstance(data, dict):
            return list(data.items())

        else:
            return [data]


class JsonWriter:
    """Write JSON files with extra features."""
    def __init__(
        self,
        filepath: PathLike,
        data: dict | list | tuple | int | float | bool | None = None,
        numpy_safe: bool = True,
        **kwargs: tp.Any
    ) -> None:
        """
        Write JSON files with extra features.

        Parameters
        ----------
        filepath: str | Path
            Path to the file to write.

        data: dict | list | tuple | int | float | bool | None, optional
            Data to write in the file. Can be added after initialization.

        numpy_safe: bool
            Whether to use a more flexible JSON encoder that auto-converts
            numpy types into python native types before serialization.
            Note that any valid custom encoder passed as keyword argument 'cls'
            will take priority over any internal one.
            Defaults to True.

        kwargs: Any
            Additional keyword arguments to pass to `json.dump()`.
        """
        check_file_format(filepath, allowed_formats="json")
        self.filepath = Path(filepath)
        self.data = data
        self.encoder = NumpyJSONEncoder if numpy_safe else json.JSONEncoder

        custom_encoder = kwargs.pop("cls", None)
        if custom_encoder is not None:
            if not issubclass(custom_encoder, json.JSONEncoder):
                raise TypeError(
                    "'cls' expected a subclass of 'json.JSONEncoder', "
                    f"got {type(custom_encoder).__name__!r}."
                )
            self.encoder = custom_encoder

        self.kwargs = kwargs

    @property
    def data(self) -> tp.Any:
        return self._data

    @data.setter
    def data(self, data: tp.Any) -> None:
        if not isinstance(data, (dict, list, tuple, int, float, bool, type(None))):
            raise TypeError(
                f"'data' of type {type(data).__qualname__!r} is not JSON-serializable."
        )
        self._data = data

    @data.deleter
    def data(self) -> None:
        self._data = None

    def write(self) -> None:
        """Write the JSON data as-is."""
        self.filepath.parent.mkdir(parents=True, exist_ok=True)
        with self.filepath.open("wt", encoding="utf-8") as fp:
            json.dump(self._data, fp, cls=self.encoder, **self.kwargs)

    def write_as_dict(self) -> None:
        """
        Write the file after casting the data into a dict explicitly.
        - If it is already a dict, it is written untouched.
        - If it is a list or tuple, it will be cast into a dict with position indices
        as integer keys.
        - In any other case, the whole data is stored in a one-element dict at integer key 0.
        """
        if isinstance(self._data, dict):
            written_data = self._data

        elif isinstance(self._data, (list, tuple)):
             written_data = {k: v for k, v in enumerate(self._data)}

        else:
            written_data = {0: self._data}

        self._data = written_data
        self.write()

    def write_as_list(self) -> None:
        """
        Write the file after casting the data into a list explicitly.
        - If it is already a list or tuple, it is written untouched
        (They are converted into JSON array the same way).
        - If it is a dict, it is cast into a list of key-value pairs tuples.
        - In any other case, all the data is stored in a one-element list.
        """
        if isinstance(self._data, (list, tuple)):
            written_data = self._data

        elif isinstance(self._data, dict):
            written_data = list(self._data.items())

        else:
            written_data = [self._data]

        self._data = written_data
        self.write()


class NumpyJSONEncoder(json.JSONEncoder):
    """Encode Numpy types like native python types."""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        return super().default(obj)