"""
Read and parse concatenated minimal poscar formatted structures in a file.
"""

import os
from os import path
import re
import typing as tp

import numpy as np


class PoscarFile:
    """
    Read and parse concatenated minimal poscar formatted structures in a file.
    Each poscar structure has to be separated by a comment line starting with "#".
    """
    filename: str
    data: list[str]
    matcher: re.Pattern = re.compile(r"^#.*?$(?=\n#|\Z)", re.MULTILINE | re.DOTALL)

    # Data caches
    _headers_cache: list[str] | None = None
    _scale_factors_cache: list[float] | None = None
    _lattices_cache: list[np.ndarray] | None = None
    _elements_cache: list[list[str]] | None = None
    _elts_counts_cache: list[list[int]] | None = None
    _positions_basis_cache: list[str] | None
    _positions_cache: list[np.ndarray] | None = None
    _cache_list: list[str] = list(
        filter(
            lambda attr: bool(re.fullmatch(r"^_[a-z_]+_cache$", attr)),
            __dict__.keys()
        )
    )

    def __init__(self, filename: str, cache: bool = True, cache_all: bool = False) -> None:
        """
        Read and parse concatenated minimal poscar formatted structures in a file.
        Each poscar structure has to be separated by a comment line starting with "#".

        Parameters
        ----------
        filename: str
            Name of the file to parse.

        cache: bool
            Whether to cache properties values after their first call.
            Defaults to True.
        
        cache_all: bool
            Whether to compute and cache all properties at initialization.
            Defaults to False.
        """
        if not path.isfile(filename):
            raise FileNotFoundError(f"{filename}: No such file found.")

        self.filename = filename
        self.cache = cache
        self.data = self._parse_file_data(filename)

        if cache_all:
            self._headers_cache = self.headers
            self._scale_factors_cache = self.scale_factors
            self._lattices_cache = self.lattices
            self._elements_cache = self.elements
            self._elts_counts_cache = self.elts_counts
            self._positions_basis_cache = self.positions_basis
            self._positions_cache = self.positions

    def __len__(self) -> int:
        return len(self.data)
    
    def _parse_file_data(self, filename: str) -> list[str]:
        """Read the file and parse poscar structures."""
        with open(filename, "rt", encoding="utf-8") as fp:
            data = fp.read()

        return self.matcher.findall(data)

    @property
    def headers(self) -> list[str]:
        """Get the header comment line of each structure."""
        if self._headers_cache is not None:
            return self._headers_cache

        parsed = [poscar.splitlines()[0] for poscar in self.data]

        if self.cache:
            self._headers_cache = parsed

        return parsed

    @property
    def scale_factors(self) -> list[float]:
        """Get the lattice scale factor of each structure."""
        if self._scale_factors_cache is not None:
            return self._scale_factors_cache
        return [float(poscar.splitlines()[1]) for poscar in self.data]
    
    @property
    def lattices(self) -> list[np.ndarray]:
        """Get the lattice matrix in a numpy 3x3 array of each structure."""
        if self._lattices_cache is not None:
            return self._lattices_cache
        raise NotImplementedError
    
    @property
    def elements(self) -> list[list[str]]:
        """Get the element list of each structure."""
        if self._elements_cache is not None:
            return self._elements_cache
        raise NotImplementedError

    @property
    def elts_counts(self) -> list[list[int]]:
        """Get the count of each element of each structure."""
        if self._elts_counts_cache is not None:
            return self._elts_counts_cache
        raise NotImplementedError
    
    @property
    def positions_basis(self) -> list[str]:
        """
        Get the atomic positions base (i.e. 'direct' or 'cartesian') of each structure.
        """
        if self._positions_basis_cache is not None:
            return self._positions_basis_cache
        raise NotImplementedError

    @property
    def positions(self) -> list[np.ndarray]:
        """Get atomic positions of all atoms in each structure."""
        if self._positions_cache is not None:
            return self._positions_cache
        raise NotImplementedError

    def as_dict(self) -> dict[str, list[tp.Any]]:
        """Get all the data in a dict format."""
        dct = {}
        for cache_attr in self._cache_list:
            attr = cache_attr.lstrip("_").replace("_cache", "")
            dct[attr] = getattr(self, attr)
        return dct
