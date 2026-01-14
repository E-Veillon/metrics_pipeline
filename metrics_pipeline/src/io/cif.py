"""I/O operations for concatenated CIF data files and strings."""

import os
import re
from pathlib import Path
import functools as ft
import typing as tp
import typing_extensions as tpe

from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser, CifWriter, CifBlock
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmOp

from .io_base import PathLike
from src.utils import VisualIterator, raise_or_warn


class CIFParsingError(Exception):
    """An error occurred while trying to parse a CIF string."""


class CIFFile:
    """
    Load, write and parse concatenated CIF strings into structures.
    """
    matcher: re.Pattern = re.compile(r"^data_.*?$(?=\ndata_|\Z)", re.MULTILINE | re.DOTALL)

    def __init__(
        self,
        data: list[str] | None = None,
        special_keys: list[str] | None = None,
        workers: int | None = None
    ) -> None:
        """
        Load, write and parse concatenated CIF strings into structures.
        Stored CIFs are purposefully not checked unless conversion to structure is needed
        to keep the representation as close as possible to what was in the loaded data.
        The constructor is typically not called directly when loading external data
        (prefer using either `from_str()` or `from_file()`), but can be called without
        specifying `data` to create an empty instance to fill manually with CIF strings and
        write a file out of them afterward.

        Parameters
        ----------
        data: list[str], optional
            List of CIF strings, each representing an individual structure.

        special_keys: list[str], optional
            Keys pointing to data inside the CIFs that should be stored inside corresponding
            structure's `properties` dict attribute, either when converting string to structure
            or when searching for stored data when converting structure back to string. Pass
            'header' to store the header line (without the 'data' part) inside a 'header' key.

        workers: int, optional
            Number of processes to use in parallel. If not given, will use default of
            `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
            and execute sequentially.
        """
        if data is None:
            data = []

        self._data = data
        self.special_keys = special_keys
        self.workers = workers

    def _get_chunksize(self) -> int:
        """
        Compute chunksize for multiprocess tasks based on current stored data length.
        """
        return min(10, len(self._data) // 100 + 1)

    def get_cifs(self) -> list[str]:
        """Get the list of individual CIF strings."""
        return self._data

    def add_cifs(self, cifs: str | tp.Sequence[str]) -> None:
        """
        Add CIF strings to the file.
        If a single string is given, it is assumed to be a single CIF string.
        If a list or tuple is given, each element is a CIF string.
        """
        match cifs:
            case str():
                self._data.append(cifs)
            case list() | tuple() if all(isinstance(cif, str) for cif in cifs):
                self._data.extend(cifs)
            case _:
                types = (
                    f"'Sequence[{' | '.join(sorted(set(type(cif).__name__ for cif in cifs)))}]'"
                    if isinstance(cifs, tp.Sequence) else f"{type(cifs).__name__!r}"
                )
                raise TypeError(
                    f"{self.__class__.__name__}: "
                    "'cifs' expected a type 'str' or 'Sequence[str]', "
                    f"got {types}."
                )

    def clear(self) -> None:
        """Remove all stored CIF data without changing other initialized parameters."""
        self._data = []

    def _parse_structure(self, cif: str) -> Structure | None:
        """Parse a Structure object from a CIF string. Return None if parsing fails."""
        parser = CifParser.from_str(cif)
        try:
            structure = parser.parse_structures(primitive=False, on_error="raise")[0]
        except ValueError:
            return None

        if not self.special_keys:
            return structure

        props = dict.fromkeys(self.special_keys)
        cif_dict: dict = parser.as_dict().popitem()[1]

        for key in self.special_keys:
            if key == "header":
                props[key] = parser._cif.data.popitem()[1].header
            else:
                props[key] = cif_dict.get(key)

        structure.properties.update(props)
        return structure

    def parse_structures(
        self, on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> tuple[list[Structure], list[int]]:
        """
        Get Structure objects from stored CIF strings.

        Parameters
        ----------
        on_error: 'raise', 'warn', 'ignore'
            What to do in case not all CIF strings could be parsed into structures.
            Defaults to 'warn'.

        Returns
        -------
        (list[Structure], list[int])
            List of structures successfully parsed from stored CIFs, in the same order
            as the internally stored CIFs, and list of indices inside the stored CIFs
            list where the data could not be parsed.
        """
        # Avoid multiprocess overhead for only one structure
        workers = 0 if len(self._data) == 1 else self.workers
        description = "Parsing CIFs into structures"

        if workers == 0:
            structures = [
                self._parse_structure(cif)
                for cif in VisualIterator(
                    self._data,
                    desc=description,
                    unit="parsed",
                    percent=True
                )
            ]
        else:
            structures = process_map(
                self._parse_structure,
                self._data,
                max_workers=self.workers,
                chunksize=self._get_chunksize(),
                desc=description
            )

        invalid_indices = [idx for idx, struct in enumerate(structures) if struct is None]
        structures = list(filter(None, structures))

        if invalid_indices:
            msg = (
                "Some CIFs could not be parsed. See their 0-based positions in the input file "
                f"hereafter: {', '.join(sorted(map(str, invalid_indices)))}."
            )
            raise_or_warn(on_error, CIFParsingError, msg)

        return structures, invalid_indices

    @staticmethod
    def _add_symmetry_data(
        block: CifBlock,
        sym_struct: SymmetrizedStructure,
        significant_figures: int = 8,
        remove_equivalents: bool = False
    ) -> CifBlock:
        """Add symmetry infos in the CIF data of a structure."""
        block.data["_symmetry_space_group_name_H-M"] = sym_struct.spacegroup.int_symbol
        block.data["_symmetry_Int_Tables_number"] = sym_struct.spacegroup.int_number
        ops: list[SymmOp] = [op.as_xyz_str() for op in sym_struct.spacegroup]
        block.data["_symmetry_equiv_pos_site_id"] = [f"{i}" for i in range(1, len(ops) + 1)]
        block.data["_symmetry_equiv_pos_as_xyz"] = ops

        if remove_equivalents:
            # Same behavior as CifWriter with defined symprec
            unique_sites = [
                (
                    min(sites, key=lambda site: tuple(abs(x) for x in site.frac_coords)),
                    len(sites),
                )
                for sites in sym_struct.equivalent_sites
            ]
            sorted_unique_sites = sorted(
                unique_sites,
                key=lambda t: (
                    t[0].species.average_electroneg,
                    -t[1],
                    t[0].a,
                    t[0].b,
                    t[0].c,
                ),
            )
            atom_site_type_symbol = []
            atom_site_symmetry_multiplicity = []
            atom_site_fract_x = []
            atom_site_fract_y = []
            atom_site_fract_z = []
            atom_site_label = []
            atom_site_occupancy = []
            for idx, (site, mult) in enumerate(sorted_unique_sites):
                for species, occupancy in site.species.items():
                    atom_site_type_symbol.append(str(species))
                    atom_site_symmetry_multiplicity.append(f"{mult}")
                    atom_site_fract_x.append(f"{site.a:.{significant_figures}f}")
                    atom_site_fract_y.append(f"{site.b:.{significant_figures}f}")
                    atom_site_fract_z.append(f"{site.c:.{significant_figures}f}")
                    atom_site_label.append(
                        f"{species.symbol}{idx}"
                        if site.label == site.species_string else site.label
                    )
                    atom_site_occupancy.append(str(occupancy))

            block.data["_atom_site_type_symbol"] = atom_site_type_symbol
            block.data["_atom_site_label"] = atom_site_label
            block.data["_atom_site_symmetry_multiplicity"] = atom_site_symmetry_multiplicity
            block.data["_atom_site_fract_x"] = atom_site_fract_x
            block.data["_atom_site_fract_y"] = atom_site_fract_y
            block.data["_atom_site_fract_z"] = atom_site_fract_z
            block.data["_atom_site_occupancy"] = atom_site_occupancy

        return block

    def _convert_to_cif(
        self,
        structure: Structure | SymmetrizedStructure,
        symmetrize: bool = False,
        symprec: float = 0.1,
        angleprec: float = 5.0,
        refine: bool = False,
        remove_equivalents: bool = False,
        significant_figures: int = 8
    ) -> str:
        """
        Convert a structure to CIF formatted string. See the `add_structures()` method
        for informations on the arguments. Do not use this method directly.
        """
        # Symmetrize and refine structure
        if symmetrize and not isinstance(structure, SymmetrizedStructure):
            spga = SpacegroupAnalyzer(structure, symprec, angleprec)

            if refine:
                structure = SpacegroupAnalyzer(
                    spga.get_refined_structure(), symprec, angleprec
                ).get_symmetrized_structure()
            else:
                structure = spga.get_symmetrized_structure()

        # Get unsymmetrized auto-generated CifBlock of the structure
        cif_block = CifWriter(
            structure, significant_figures=significant_figures
        ).cif_file.data.popitem()[1]

        # Add eventual symmetry informations manually
        if isinstance(structure, SymmetrizedStructure):
            cif_block = self._add_symmetry_data(
                cif_block, structure, significant_figures, remove_equivalents
            )
        
        # Add back saved special data, before actual structure description
        if self.special_keys:
            saved_data = []
            for key in self.special_keys:
                if key == "header":
                    cif_block.header = structure.properties[key]
                else:
                    saved_data.append((key, structure.properties[key]))
            saved_data.extend(list(cif_block.data.items()))
            cif_block.data = dict(saved_data)
        
        return str(cif_block)

    def add_structures(
        self, structures: list[Structure | SymmetrizedStructure], *args, **kwargs
    ) -> None:
        """
        Add structures data into the file, with multiprocess if `workers` was not set
        to 0 at initialization.

        Parameters
        ----------
        structures: list[Structure | SymmetrizedStructure]
            List of structures to add to the file.

        symmetrize: bool
            Whether to search for the spacegroup symmetry before converting to string.
            If a SymmetrizedStructure object is passed, this step is ignored for efficiency.
            Defaults to False.

        symprec: float
            Tolerance for symmetry finding. See pymatgen's `SpacegroupAnalyzer` for more details.
            Defaults to 0.1 (instead of pymatgen's 0.01) for better computation stability with
            simulated crystals.

        angleprec: float
            Angle tolerance for symmetry finding. Defaults to 5.0 degrees.

        refine: bool
            Whether sites inside the newly symmetrized structure should be moved to their
            expected positions according to spacegroup. Ignored if `symmetrize` is `False`
            or if a SymmetrizedStructure is passed. Defaults to False.

        remove_equivalents: bool
            Whether to only keep symmetrically distinct sites in the final string.
            **WARNING**: the final string may generate a different structure if converted again
            with the `parse_structures()` method, as it assumes all sites are present in data.
            Defaults to False.

        significant_figures: int
            Number of decimals to keep for floats formatting. Defaults to 8.
        """
        # Avoid multiprocess overhead for only one structure
        workers = 0 if len(structures) == 1 else self.workers
        description = "Converting structures to CIF"

        if workers == 0:
                cif_strings = [
                    self._convert_to_cif(structure, *args, **kwargs)
                    for structure in VisualIterator(
                        structures,
                        desc=description,
                        unit="converted",
                        percent=True
                    )
                ]
        else:
            convert_to_cif = ft.partial(self._convert_to_cif, *args, **kwargs)
            cif_strings = process_map(
                convert_to_cif,
                structures,
                max_workers=workers,
                chunksize=self._get_chunksize(),
                desc=description
            )
        self._data.extend(cif_strings)

    @classmethod
    def from_str(cls, string: str, **kwargs) -> tpe.Self:
        """
        Load CIFs from a concatenated string.

        Parameters
        ----------
        string: str
            String to parse.

        kwargs: Any
            Additional keyword arguments to pass to the constructor.
        """
        return cls(
            [str(cif).strip() for cif in re.findall(cls.matcher, string)],
            **kwargs
        )

    @classmethod
    def from_file(cls, filename: PathLike, **kwargs) -> tpe.Self:
        """
        Load a concatenated CIF file.
        
        Parameters
        ----------
        filename: str | Path
            File to load.

        kwargs: Any
            Additional keyword arguments to pass to the constructor.
        """
        with open(filename, "rt", encoding="utf-8") as fp:
            str_data = fp.read()

        return cls.from_str(str_data, **kwargs)

    def __getitem__(self, index: int) -> str:
        """Get data for one structure."""
        return self._data[index]

    def __str__(self) -> str:
        """File content as a single string."""
        return "\n".join(self._data)

    def __len__(self) -> int:
        """Number of stored CIF strings."""
        return len(self._data)

    def __eq__(self, other: tpe.Self) -> bool:
        """Whether two files contain exactly the same data."""
        if not isinstance(other, CIFFile):
            return NotImplemented

        return self is other or str(self) == str(other)

    def write_file(self, filename: PathLike) -> None:
        """Write concatenated CIF strings into a file."""
        os.makedirs(Path(filename).parent, exist_ok=True)

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write(str(self))
