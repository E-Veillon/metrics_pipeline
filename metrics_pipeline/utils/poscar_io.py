"""
Read and parse concatenated minimal poscar formatted structures in a file.
"""

import os
from os import path
import re
import typing as tp
from dataclasses import dataclass

import numpy as np


@dataclass
class PoscarBlock:
    header: str
    scale_factor: float
    lattice: np.ndarray
    elements: list[str]
    elts_count: list[int]
    positions_basis: str
    positions: np.ndarray


class PoscarFile:
    """
    Read and parse concatenated minimal VASP 5.0+ poscar formatted structures in a file.
    Minimal means each data block only contains a comment header, lattice, elements and
    positions informations.
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
    _positions_basis_cache: list[str] | None = None
    _positions_cache: list[np.ndarray] | None = None

    def __init__(self, filename: str, cache: bool = True, cache_all: bool = False) -> None:
        """
        Read and parse concatenated minimal VASP 5.0+ poscar formatted structures in a file.
        Minimal means each data block only contains a comment header, lattice, elements and
        positions informations.
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

        self._cache_list: list[str] = list(
            filter(
                lambda attr: bool(re.fullmatch(r"\A_[a-z_]+_cache\Z", attr)),
                self.__dict__.keys()
            )
        )

        self.filename = filename
        self.cache = cache
        self.cache_all = cache_all
        self.data = self._parse_file_data(filename)

        if self.cache_all:
            for cache_attr in self._cache_list:
                attr_value = getattr(self, self._get_attr_from_cache(cache_attr))
                setattr(self, cache_attr, attr_value)

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, idx: int) -> PoscarBlock:
        return PoscarBlock(
            header=self.headers[idx],
            scale_factor=self.scale_factors[idx],
            lattice=self.lattices[idx],
            elements=self.elements[idx],
            elts_count=self.elts_counts[idx],
            positions_basis=self.positions_basis[idx],
            positions=self.positions[idx]
        )

    def _parse_file_data(self, filename: str) -> list[str]:
        """Read the file and parse poscar structures."""
        with open(filename, "rt", encoding="utf-8") as fp:
            data = fp.read()

        return self.matcher.findall(data)

    def _get_attr_from_cache(self, cache_attr: str) -> str:
        return cache_attr.lstrip("_").replace("_cache", "")

    @property
    def headers(self) -> list[str]:
        """Get the header comment line of each structure."""
        if self._headers_cache is not None:
            return self._headers_cache

        parsed = [poscar.splitlines()[0][1:].strip() for poscar in self.data]

        if self.cache and not self.cache_all:
            self._headers_cache = parsed

        return parsed

    @property
    def scale_factors(self) -> list[float]:
        """Get the lattice scale factor of each structure."""
        if self._scale_factors_cache is not None:
            return self._scale_factors_cache

        parsed = [float(poscar.splitlines()[1]) for poscar in self.data]

        if self.cache and not self.cache_all:
            self._scale_factors_cache = parsed

        return parsed

    @property
    def lattices(self) -> list[np.ndarray]:
        """Get the lattice matrix in a numpy 3x3 array of each structure."""
        if self._lattices_cache is not None:
            return self._lattices_cache

        parsed = [
            np.array(
                [
                    list(map(float, poscar.splitlines()[2])),
                    list(map(float, poscar.splitlines()[3])),
                    list(map(float, poscar.splitlines()[4]))
                ], dtype=np.float32
            ) for poscar in self.data
        ]

        if self.cache and not self.cache_all:
            self._lattices_cache = parsed

        return parsed
    
    @property
    def elements(self) -> list[list[str]]:
        """Get the element list of each structure."""
        if self._elements_cache is not None:
            return self._elements_cache

        parsed = [poscar.splitlines()[5].split() for poscar in self.data]

        if self.cache and not self.cache_all:
            self._elements_cache = parsed

        return parsed

    @property
    def elts_counts(self) -> list[list[int]]:
        """Get the count of each element of each structure."""
        if self._elts_counts_cache is not None:
            return self._elts_counts_cache

        parsed = [list(map(int, poscar.splitlines()[6].split())) for poscar in self.data]

        if self.cache and not self.cache_all:
            self._elts_counts_cache = parsed

        return parsed
    
    @property
    def positions_basis(self) -> list[str]:
        """
        Get the atomic positions base (i.e. 'direct' or 'cartesian') of each structure.
        """
        if self._positions_basis_cache is not None:
            return self._positions_basis_cache

        parsed = [poscar.splitlines()[7] for poscar in self.data]

        if self.cache and not self.cache_all:
            self._positions_basis_cache = parsed

        return parsed

    @property
    def positions(self) -> list[np.ndarray]:
        """Get atomic positions of all atoms in each structure."""
        if self._positions_cache is not None:
            return self._positions_cache

        parsed = []
        for poscar in self.data:
            positions = list(map(float, poscar.splitlines()[8:]))
            assert len(positions) % 3 == 0, (
                "The number of position values is not a multiple of 3 in the structure below:\n"
                f"{poscar}\n"
            )
            atoms_pos = []
            while positions:
                atoms_pos.append(positions[0:3])
                positions = positions[3:]
            parsed.append(np.array(atoms_pos, dtype=np.float32))

        if self.cache and not self.cache_all:
            self._positions_cache = parsed

        return parsed

    def as_dict(self) -> dict[str, list[tp.Any]]:
        """Get all the data in a dict format."""
        dct = {}
        for cache_attr in self._cache_list:
            attr = self._get_attr_from_cache(cache_attr)
            dct[attr] = getattr(self, attr)
        return dct
