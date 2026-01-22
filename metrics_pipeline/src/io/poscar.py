"""
Read and parse concatenated minimal POSCAR formatted structures from a file.
"""

import os
import re
import typing as tp
import typing_extensions as tpe
from dataclasses import dataclass, asdict
import itertools as itt

import numpy as np
import numpy.typing as npt

from .io_base import PathLike
from src.utils import ALL_ELT_SYMBOL_TO_Z


@dataclass
class PoscarBlock:
    """
    Convenient dataclass for storing minimal VASP 5.0+ POSCAR data of a structure.

    Attributes
    ----------
    header: str
        Comment line at the beginning.

    scale_factor: float
        Global scaling factor applied to lattice vectors and cartesian atomic positions.

    lattice: ndarray
        A 3x3 numpy array of the lattice vectors components in cartesian space.

    elements: list[str]
        List of chemical elements symbols.

    elts_count: list[int]
        List of the number of atoms in the structure of each element.
        Must be in the same order as 'elements' for correct mapping.

    positions_basis: str
        Basis of atomic positions, either 'direct' for fractional coordinates or 'cartesian'
        for cartesian coordinates.

    positions: ndarray
        Nx3 array of atomic positions, N being the number of atoms.
    """
    header: str
    scale_factor: float
    lattice: npt.NDArray[np.float32]
    elements: list[str]
    elts_count: list[int]
    positions_basis: str
    positions: npt.NDArray[np.float32]

    def __len__(self) -> int:
        """Number of atoms in data."""
        return self.positions.shape[0]

    def __str__(self) -> str:
        """Convenient default str formatting."""
        return self.as_string()

    @staticmethod
    def _check_string(string: str) -> None:
        """Assertions to verify if given string conforms to POSCAR formatting."""
        lines = string.strip().splitlines()
        assert len(lines) >= 9, ValueError(
            "Given string must contain at least 9 lines to be a valid POSCAR: "
            "1 header, 1 scale factor, 3 lattice vectors, 1 element list, "
            "1 element counts list, 1 position basis ('direct' or 'cartesian'), "
            "and at least 1 atomic position."
        )
        assert lines[0].startswith("#"), ValueError(
            f"1st line must be a header comment line starting with a '#', got {lines[0]!r}."
        )
        assert is_float_string(lines[1]), ValueError(
            f"2nd line must be a floating point number, got {lines[1]!r}."
        )
        for idx, line in enumerate(lines[2:5], start=1):
            coeffs = line.split()
            all_floats = all(is_float_string(coeff) for coeff in coeffs)
            assert len(coeffs) == 3 and all_floats, ValueError(
                f"Lattice vector line {idx} must be 3 floating point numbers, got {line!r}."
            )
        valid_elements = set(ALL_ELT_SYMBOL_TO_Z)
        assert all(elt in valid_elements for elt in lines[5].split()), ValueError(
            "6th line contains data that is not a valid element symbol. "
            f"Line got: {lines[5]!r}."
        )
        assert all(s.isdecimal() for s in lines[6].split()), ValueError(
            "7th line contains data that is not an integer. "
            f"Line got: {lines[6]!r}."
        )
        assert len(lines[5].split()) == len(lines[6].split()), ValueError(
            "6th and 7th lines should have the same number of data, "
            f"got {len(lines[5].split())} for 6th line and {len(lines[6].split())} for 7th line."
        )
        assert lines[7] in {"direct", "cartesian"}, ValueError(
            f"8th line should be either 'direct' or 'cartesian', got {lines[7]!r}."
        )
        elts_list = get_element_list(lines[5].split(), list(map(int, lines[6].split())))
        assert len(elts_list) == len(lines[8:]), ValueError(
            f"Declared number of atoms ({len(elts_list)}) does not match with the number "
            f"of atomic positions ({len(lines[8:])})."
        )
        for idx, atom_line in enumerate(lines[8:]):
            position_str = atom_line.split()
            if len(position_str) == 4 and position_str[3] in valid_elements:
                assert position_str[3] == elts_list[idx], ValueError(
                    f"Element shown next to atomic position line {idx + 1} does not match "
                    "with global element order."
                )
                position_str = position_str[:3]
            assert len(position_str) == 3, ValueError(
                f"The number of position values is not 3 for atomic position line {idx + 1}."
            )

    def is_valid(self) -> bool:
        """Whether the object misses any data or contains out-of-specs data."""
        try:
            self._check_string(str(self))
        except (ValueError, AssertionError):
            return False
        return True

    # ===== I/O methods =====
    def as_string(self, decimals: int = 8):
        """
        Get POSCAR formatted string from data.

        Parameters
        ----------
        decimals: int
            Number of decimal places for lattice vectors and atomic positions.
            Defaults to 8.

        Returns
        -------
        str
            POSCAR formatted string of data.
        """
        assert decimals > 0, ValueError(
            f"{self.as_string.__qualname__}: 'decimals' argument must be > 0."
        )
        lines = [f"# {self.header}"]
        lines.append(f"{self.scale_factor}")
        BASE_RJUST = 5 # Enough space to go from -999.* to 9999.*
        RJUST = BASE_RJUST + decimals
        for lattice_line_idx in range(self.lattice.shape[0]):
            coeffs = [
                f"{coeff:>{RJUST}.{decimals}f}" for coeff in self.lattice[lattice_line_idx, :]
            ]
            lines.append(f"{' '.join(coeffs)}")
        lines.append(f"{' '.join(self.elements)}")
        lines.append(f"{' '.join(list(map(str, self.elts_count)))}")
        lines.append(f"{self.positions_basis}")
        elt_list = get_element_list(self.elements, self.elts_count)
        for pos_line_idx, elt in zip(range(self.positions.shape[0]), elt_list):
            coeffs = [
                f"{coeff:>{RJUST}.{decimals}f}" for coeff in self.positions[pos_line_idx, :]
            ]
            lines.append(f"{' '.join(coeffs)} {elt}")

        return "\n".join(lines)

    @classmethod
    def from_string(cls, string: str) -> tpe.Self:
        """Get a PoscarBlock from POSCAR-like formatted string."""
        cls._check_string(string)
        lines = string.strip().splitlines()
        header = lines[0][1:].strip()
        scale_factor = float(lines[1])
        lattice = np.array(
            [
                list(map(float, lines[2])),
                list(map(float, lines[3])),
                list(map(float, lines[4]))
            ], dtype=np.float32
        )
        elements = lines[5].split()
        elts_count = list(map(int, lines[6].split()))
        pos_basis = lines[7]
        pos_list = [list(map(float, atom_line.split())) for atom_line in lines[8:]]
        positions = np.array(pos_list, dtype=np.float32)
        return cls(header, scale_factor, lattice, elements, elts_count, pos_basis, positions)

    def as_dict(self) -> dict[str, tp.Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, dct: dict[str, tp.Any]) -> tpe.Self:
        return cls(
            header=dct["header"],
            scale_factor=dct["scale_factor"],
            lattice=dct["lattice"],
            elements=dct["elements"],
            elts_count=dct["elts_count"],
            positions_basis=dct["positions_basis"],
            positions=dct["positions"],
        )


class PoscarFile:
    """
    Read and parse concatenated minimal VASP 5.0+ poscar formatted structures in a file.
    Minimal means each data block only contains a comment header, lattice, elements and
    positions informations.
    Each poscar structure has to begin by a comment line starting with "#".
    """
    filename: PathLike | None
    _data: list[str]
    matcher: re.Pattern = re.compile(r"^#.*?$(?=\n#|\Z)", re.MULTILINE | re.DOTALL)

    # Data caches
    _headers_cache: list[str] | None = None
    _scale_factors_cache: list[float] | None = None
    _lattices_cache: list[np.ndarray] | None = None
    _elements_cache: list[list[str]] | None = None
    _elts_counts_cache: list[list[int]] | None = None
    _positions_basis_cache: list[str] | None = None
    _positions_cache: list[np.ndarray] | None = None

    def __init__(
        self,
        filename: PathLike | None = None,
        strict: bool = True,
        cache: bool = True,
        cache_all: bool = False
    ) -> None:
        """
        Read and parse concatenated minimal VASP 5.0+ poscar formatted structures in a file.
        Minimal means each data block only contains a comment header, lattice, elements and
        positions informations.
        Each poscar structure has to be separated by a comment line starting with "#".

        Parameters
        ----------
        filename: str, optional
            Name of the file to parse. If not given, creates an empty instance that can be
            manually populated.

        strict: bool
            Whether to verify validity of each structure during data parsing.
            Defaults to True.

        cache: bool
            Whether to cache properties values after their first call.
            Defaults to True.

        cache_all: bool
            Whether to compute and cache all properties at initialization.
            Defaults to False.
        """
        self._cache_list: list[str] = list(
            filter(
                lambda attr: bool(re.fullmatch(r"\A_[a-z_]+_cache\Z", attr)),
                self.__dict__.keys()
            )
        )
        if filename is None:
            self._create_empty_instance()

        else:
            if not os.path.isfile(filename):
                    raise FileNotFoundError(f"{filename}: No such file found.")

            self.filename = filename
            self.cache = cache
            self.cache_all = cache_all
            self._data = self._parse_file_data(filename, strict)

            if self.cache_all:
                for cache_attr in self._cache_list:
                    attr_value = getattr(self, self._get_attr_from_cache(cache_attr))
                    setattr(self, cache_attr, attr_value)

    def __len__(self) -> int:
        return len(self._data)

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

    # ===== Internal helper methods =====
    @staticmethod
    def _check_data(data: list[str]) -> None:
        """Verify validity of each parsed structure string."""
        for idx, block in enumerate(data):
            try:
                PoscarBlock._check_string(block)
            except (ValueError, AssertionError) as exc:
                raise ValueError(
                    f"An error occurred while trying to parse structure {idx}. "
                    f"See below for details:\n{exc}"
                )

    def _parse_file_data(self, filename: PathLike, strict: bool = True) -> list[str]:
        """Read the file and parse poscar structures."""
        with open(filename, "rt", encoding="utf-8") as fp:
            file_data = fp.read()

        data = [str(block).strip() for block in self.matcher.findall(file_data)]

        if strict:
            self._check_data(data)

        return data

    def _get_attr_from_cache(self, cache_attr: str) -> str:
        return cache_attr.lstrip("_").replace("_cache", "")

    def _create_empty_instance(self) -> None:
        """Create empty instance to populate manually."""
        self.filename = None
        self.cache = True
        self.cache_all = True
        self._data = []

    # ===== Public methods and properties =====
    def is_valid(self) -> bool:
        """Whether the object misses any data or contains out-of-specs data."""
        return all(block.is_valid() for block in self)

    def add_poscar_block(self, block: PoscarBlock, decimals: int = 8) -> None:
        """
        Append the data from a PoscarBlock object at the end of the data of this object in place.
        Adding valid PoscarBlock to a PoscarFile is the safest way to avoid missing data.

        Parameters
        ----------
        block: PoscarBlock
            PorcarBlock object to add to the data.

        decimals: int
            Number of decimals to keep in the string representation of the PoscarBlock
            for lattice vectors and atomic positions. Defaults to 8.
        """
        assert block.is_valid(), ValueError(
            "Given PoscarBlock is not filled with valid data only. "
            "It is not added to the PoscarFile in order to avoid file's data corruption."
        )
        self._data.append(block.as_string(decimals))
        for cache_attr in self._cache_list:
            attr_content: list[tp.Any] | None = getattr(self, cache_attr)
            if attr_content is None:
                continue
            attr = self._get_attr_from_cache(cache_attr)
            # Remove plural 's' for some attributes to correspond
            block_attr = attr if hasattr(block, attr) else attr[:-1]
            assert hasattr(block, block_attr), AttributeError(
                f"Attribute {block_attr}(s) in {self.__class__.__name__} does not match "
                f"any {block.__class__.__name__} attribute. This error is likely caused "
                "by a bug in the code and should be reported as issue to the developpers."
            )
            attr_content.append(getattr(block, block_attr))
            setattr(self, cache_attr, attr_content)

    @property
    def data(self) -> list[str]:
        """Get the POSCAR formatted string of each structure."""
        return self._data

    @property
    def headers(self) -> list[str]:
        """Get the header comment line of each structure."""
        if self._headers_cache is not None:
            return self._headers_cache

        parsed = [poscar.splitlines()[0][1:].strip() for poscar in self._data]

        if self.cache and not self.cache_all:
            self._headers_cache = parsed

        return parsed

    @property
    def scale_factors(self) -> list[float]:
        """Get the lattice scale factor of each structure."""
        if self._scale_factors_cache is not None:
            return self._scale_factors_cache

        parsed = [float(poscar.splitlines()[1]) for poscar in self._data]

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
                    list(map(float, poscar.splitlines()[2].split())),
                    list(map(float, poscar.splitlines()[3].split())),
                    list(map(float, poscar.splitlines()[4].split()))
                ], dtype=np.float32
            ) for poscar in self._data
        ]

        if self.cache and not self.cache_all:
            self._lattices_cache = parsed

        return parsed

    @property
    def elements(self) -> list[list[str]]:
        """Get the element list of each structure."""
        if self._elements_cache is not None:
            return self._elements_cache

        parsed = [poscar.splitlines()[5].split() for poscar in self._data]

        if self.cache and not self.cache_all:
            self._elements_cache = parsed

        return parsed

    @property
    def elts_counts(self) -> list[list[int]]:
        """Get the count of each element of each structure."""
        if self._elts_counts_cache is not None:
            return self._elts_counts_cache

        parsed = [list(map(int, poscar.splitlines()[6].split())) for poscar in self._data]

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

        parsed = [poscar.splitlines()[7] for poscar in self._data]

        if self.cache and not self.cache_all:
            self._positions_basis_cache = parsed

        return parsed

    @property
    def positions(self) -> list[np.ndarray]:
        """Get atomic positions of all atoms in each structure."""
        if self._positions_cache is not None:
            return self._positions_cache

        parsed = []
        for poscar in self._data:
            lines = poscar.splitlines()
            atoms_pos = [list(map(float, atom_line.split()[:3])) for atom_line in lines[8:]]
            parsed.append(np.array(atoms_pos, dtype=np.float32))

        if self.cache and not self.cache_all:
            self._positions_cache = parsed

        return parsed

    # ===== I/O methods =====
    def as_dict(self) -> dict[str, list[tp.Any]]:
        """Get all the data in a dict format."""
        dct = {}
        for cache_attr in self._cache_list:
            attr = self._get_attr_from_cache(cache_attr)
            dct[attr] = getattr(self, attr)
        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, list[tp.Any]], strict: bool = True) -> tpe.Self:
        """
        Regenerate a PoscarFile object from a dict of same format
        as given by the `as_dict()` method.

        Parameters
        ----------
        dct: dcit[str, list]
            Dict of lists of data. Keys must map exactly to attribute names.

        strict: bool
            Whether to verify data validity while parsing. Defaults to True.

        Returns
        -------
        PoscarFile
            Parsed PoscarFile object.
        """
        pfile = cls()
        for cache_attr in pfile._cache_list:
            attr = pfile._get_attr_from_cache(cache_attr)
            assert dct.get(attr, False), KeyError(
                f"{cls.from_dict.__qualname__}: dictionary lacks"
                f"the following key or it is empty: {attr!r}."
            )
            setattr(pfile, cache_attr, dct[attr])

        data = [block.as_string() for block in pfile]

        if strict:
            pfile._check_data(data)

        pfile._data = data

        return pfile

    def write_file(self, filename: str, exist_ok: bool = False) -> None:
        """
        Write formatted data into a file.
        All data must be present and follow specs for all structures for this to work.

        Parameters
        ----------
        filename: str
            Path to write the file.

        exist_ok: bool
            If given path already exists, whether to override it or raise a `FileExistsError`.
            Defaults to False, which will raise an error.

        Raises
        ------
        FileExistsError
            If given path already exists and `exist_ok` is set to False.

        MissingDataError
            If some data in the instance is missing or not valid.
        """
        assert self.is_valid(), ValueError(
            "This instance contains some missing or out-of-specs data, "
            "preventing proper file writing."
        )
        os.makedirs(os.path.dirname(filename), exist_ok=exist_ok)
        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(self._data))


def is_float_string(string: str) -> bool:
    """Check if given string is a floating point number."""
    parts = string.split(".")
    return len(parts) == 2 and all(s.lstrip("-").isdecimal() for s in parts)


def get_element_list(elements: list[str], counts: list[int]) -> list[str]:
    """Get a list of elements symbols in order of appearance."""
    return list(itt.chain.from_iterable(
            [[elt] * num_elt for elt, num_elt in zip(elements, counts)]
        ))
