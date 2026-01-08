"""Management of VASP files for quantum simulations."""

import os
import re
from pathlib import Path
import typing as tp
from collections import defaultdict
from collections.abc import Callable
import json
import warnings
import itertools as itt
import functools as ft
from enum import Enum

from monty.os import cd
import xml.etree.ElementTree as ET
from tqdm.contrib.concurrent import process_map

from pymatgen.core import Structure
from pymatgen.io.vasp import VaspInput, Vasprun, Poscar, Xdatcar

from .io_base import PathLike, check_file_or_dir
from src.utils import raise_or_warn, VisualIterator


class VaspParsingError(Exception):
    """An error occurred when trying to parse a VASP output file."""


class VaspWriter:
    """
    Lightweight wrapper with strict assertions to write multiple VASP run directories at once.
    """
    def __init__(
        self,
        base_dir: PathLike,
        vasp_inputs: dict[str, VaspInput] | None = None,
    ) -> None:
        """
        Lightweight wrapper to write multiple VASP run directories at once.

        Parameters
        ----------
        base_dir: PathLike
            Path to base directory to write VASP run directories to. It is recommended to
            use an empty directory. Creates the path if it does not exist.

        vasp_inputs: dict[str, VaspInput], optional
            Dict of {str: VaspInput} with VASP run directories names to create as keys and
            VaspInput objects containing all necessary files to write to the run directories
            as values. They can be added directly at initialization and / or one by one later
            with the `add_run_dir()` method. An extra safety check in the method prevents
            overriding an already registered run by providing the same name twice.
        """
        os.makedirs(base_dir, exist_ok=True)
        self.base_dir = Path(base_dir)
        
        if vasp_inputs is None:
            self.vasp_inputs = {}

        else:
            assert isinstance(self.vasp_inputs, dict), TypeError(
                f"'vasp_inputs' expected a type 'dict', got {type(self.vasp_inputs).__name__}."
            )
            for run_name, vasp_input in vasp_inputs.items():
                assert isinstance(run_name, str), TypeError(
                    f"'vasp_inputs' key {run_name} must be of type 'str', "
                    f"got {type(run_name).__name__}."
                )
                assert isinstance(vasp_input, VaspInput), TypeError(
                    f"'vasp_inputs[{run_name}]' value must be of VaspInput type, "
                    f"got {type(vasp_input).__name__}."
                )
                self._check_vasp_input(run_name, vasp_input)
        
            self.vasp_inputs = vasp_inputs

    @staticmethod
    def _check_vasp_input(run_name: str, vasp_input: VaspInput) -> None:
        """Assert that all mandatory files are present in the VaspInput."""
        assert vasp_input.get("INCAR") is not None, ValueError(
                f"{run_name}: There is no INCAR defined in the input !"
        )
        assert vasp_input.get("POSCAR") is not None, ValueError(
                f"{run_name}: There is no POSCAR defined in the input !"
        )
        assert (
            vasp_input.get("KPOINTS") is not None or
            vasp_input["INCAR"].get("KSPACING") is not None
        ), ValueError(
            f"{run_name}: There is no KPOINTS or KSPACING tag defined in the input !"
        )
        assert vasp_input.get("POTCAR") is not None, ValueError(
                f"{run_name}: There is no POTCAR defined in the input !"
        )

    def add_run_dir(
        self, run_name: str, vasp_input: VaspInput, check_override: bool = True
    ) -> None:
        """
        Add a run directory to the writer.
        
        Parameters
        ----------
        run_name: str
            Name of the new run directory.

        vasp_input: VaspInput
            VaspInput object containing data to write VASP input files in the new run directory.

        check_override: bool
            Whether to check if run name is already registered in the writer before overriding it.
            Defaults to True.
        """
        self._check_vasp_input(run_name, vasp_input)
        if check_override:
            assert self.vasp_inputs.get(run_name) is None, ValueError(
                f"{run_name}: This run already exists in registered runs."
            )
        self.vasp_inputs[run_name] = vasp_input

    def write_run_dirs(self) -> None:
        """
        Write all initialized run directories and input files.
        """
        for run_name, vasp_input in self.vasp_inputs.items():
            run_dir = os.path.join(self.base_dir, run_name)
            vasp_input.write_input(output_dir=run_dir)


class VaspParser:
    """Match and parse VASP run directories."""
    def __init__(
        self,
        base_dir: PathLike,
        indices: list[int] | None = None,
        max_index: int | None = None,
        match_all: bool = True,
        unique: bool = True,
        workers: int | None = None
    ) -> None:
        """
        Match and parse VASP run directories.
        
        Parameters
        ----------
        base_dir: PathLike
            Base directory containing VASP finished runs to extract or to write
            new VASP run directories to.

        indices: list[int], optional
            Indices of wanted structure directories to match.
            If not given, defaults to all possible indices between 0 and `max_index`.

        max_index: int, optional
            Maximum index to try to match if `indices` is not given. If not given,
            defaults to the max found index in the screened directory.

        match_all: bool
            Whether each index should match at least one directory.
            If True, raise an error when an index does not match any directory.
            Defaults to True.

        unique: bool
            Whether each index should match at most one directory.
            If True, raise an error when an index matches several directories.
            Defaults to True.

        workers: int, optional
            Number of processes to use in parallel. If not given, will use default of
            `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
            and execute sequentially.

        Raises
        ------
        - `FileNotFoundError` if match_all is True and an index does not match any directory.
        - `FileExistsError` if unique is True and an index matches several directories.

        Notes
        -----
        - The default behaviour of `indices` only goes up to the max found index, meaning that
        if higher indices should have been present they will not be caught by the `match_all`
        check. If properly matching all indices in a known range is what you need, provide
        wanted max to `max_index`.
        """
        check_file_or_dir(base_dir, "dir", check_empty=True)
        self.base_dir = Path(base_dir)

        # Match all existing structure directories in base directory
        subtree = list(self.base_dir.iterdir())
        matching_dirs = sorted(filter(self.is_struct_dir, subtree))
        struct_indices = [int(dir.name.split("_")[0]) for dir in matching_dirs]

        struct_dirs_dict = defaultdict(list)
        for idx, dir in zip(struct_indices, matching_dirs):
            struct_dirs_dict[idx].append(dir)

        # Setup indices
        if indices is None:
            indices = list(
                range(max(struct_indices) + 1) if max_index is None else range(max_index + 1)
            )
        else:
            assert all(isinstance(idx, int) for idx in indices)

        self.indices = sorted(indices)
        self.match_all = match_all
        self.unique = unique
        self.workers = workers
        self._all_struct_dirs = struct_dirs_dict
        self._struct_dirs = self._match_struct_dirs_indices()

    @staticmethod
    def is_struct_dir(path: PathLike, strict: bool = True) -> bool:
        """
        Checks whether given path is a directory having naming convention of a structure directory.

        Parameters
        ----------
        path: str | Path
            File to check.

        strict: bool
            Use a stricter matching pattern that does not authorize letters that never appear
            in the symbols of the 118 elements (i.e. 'J', 'Q', 'j', 'q', 'w', 'x', 'z').
            Defaults to True.

        Returns
        -------
        bool
        - `True` if given path is a directory following structure directory naming.
        - `False` otherwise.
        """
        # NOTE: Here 'path' must be a Path object, do not change for os.path.
        path = Path(path)
        elt_w_idx = r"[A-IK-PR-Z][a-ik-pr-vy]?\d*" if strict else r"[A-Z][a-z]?\d*"
        pattern = fr"\A\d+_(?:{elt_w_idx}|\((?:{elt_w_idx})+\)\d*)+\Z"
        return (
            path.is_dir()
            and re.fullmatch(pattern, path.name) is not None
        )

    @property
    def all_struct_dirs(self) -> list[Path]:
        """Sorted list of all structure directories in the base directory."""
        return sorted(itt.chain.from_iterable(self._all_struct_dirs.values()))

    @property
    def struct_dirs(self) -> list[Path]:
        """Sorted list of structure directories matching initialized indices."""
        return sorted(self._struct_dirs)

    def _match_struct_dirs_indices(self) -> list[Path]:
        """Match structure directories corresponding to initialized indices."""
        if self.match_all:
            no_match_indices = [
                idx for idx in self.indices if len(self._all_struct_dirs.get(idx, [])) == 0
            ]
            if no_match_indices:
                raise FileNotFoundError(
                    "Following indices did not match any structure directory "
                    f"in {self.base_dir}: {no_match_indices}."
                )

        if self.unique:
            multi_match_indices = [
                idx for idx in self.indices if len(self._all_struct_dirs.get(idx, [])) > 1
            ]
            if multi_match_indices:
                raise FileExistsError(
                    "Following indices matched more than one structure directory "
                    f"in {self.base_dir}: {multi_match_indices}."
                )

        matching_dirs = list(itt.chain.from_iterable(
            [self._all_struct_dirs.get(idx, []) for idx in self.indices]
        ))
        return matching_dirs

    @staticmethod
    def get_safe_vasprun(
        run_dir: PathLike, filename: str = "vasprun.xml", converged: bool = True, **kwargs
    ) -> Vasprun | None:
        """
        Read and parse vasprun.xml file at given VASP run directory.
        Check the file integrity and whether the VASP run terminated normally.

        Parameters
        ----------
        run_dir: str | Path
            Path to the directory containing the VASP calculation.

        filename: str, optional
            Name of the .xml file containing the run data. Defaults to 'vasprun.xml',
            which is the automatic name given by VASP 5.0+ when the file is created.

        converged: bool
            Whether to only return the Vasprun object if the run is converged.
            Defaults to True.

        kwargs: Any
            Keyword arguments to pass to the Vasprun class.

        Returns
        -------
        Vasprun | None
            Vasprun object if the run is valid, None otherwise.
        """
        check_file_or_dir(run_dir, "dir")
        vasprun_path = os.path.join(run_dir, filename)
        check_file_or_dir(vasprun_path, "file", allowed_formats="xml")

        try:
            vasprun = Vasprun(filename=vasprun_path, **kwargs)
        except ET.ParseError:
            return None
        except UnicodeDecodeError:
            warnings.warn(
                f"{filename} file at {run_dir} contains "
                "unreadable characters for 'utf-8' codec.\n"
                "Associated data is therefore considered erroneous "
                "and is not parsed further."
            )
            return None

        if converged and not vasprun.converged:
            return None

        return vasprun

    def get_struct_from_run_dir(
        self,
        run_dir: PathLike,
        try_vasprun: bool = True,
        try_contcar: bool = True,
        try_xdatcar: bool = False,
        try_poscar: bool = False,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> Structure | None:
        """
        Attempt to extract the best structure possible from a VASP run directory.

        Algorithm
        ---------
        Attempt to extract a structure from a VASP run in following order,
        passing to the next step if the current one fails:
        - 1 - Tries to get the final structure from vasprun.xml
        - 2 - Tries to get the final structure from CONTCAR
        - 3 - Tries to get the last structure from XDATCAR
        - 4 - Tries to get the initial structure from POSCAR

        Parameters
        ----------
        path: str | Path
            Path to the VASP run directory.

        try_vasprun: bool
            Whether to try to get the final structure from vasprun.xml. Defaults to True.

        try_contcar: bool
            Whether to try to get the final structure from CONTCAR. Defaults to True.

        try_xdatcar: bool
            Whether to try to get the last structure from XDATCAR. Defaults to False.

        try_poscar: bool
            Whether to try to get the initial structure from POSCAR. Defaults to False.

        on_error: 'raise', 'warn', 'ignore'
            What to do when no structure could be extracted from attempted files.
            Defaults to "warn".

        Returns
        -------
        Structure | None
            A Structure object if it could be extracted, `None` if it could not and
            `on_error` is not set to "raise".

        Raises
        ------
        `VaspParsingError` if none of the attempts were successful and `on_error`
        is set to "raise".
        """
        check_file_or_dir(run_dir, "dir")

        if try_vasprun:
            try:
                file = os.path.join(run_dir, "vasprun.xml")
                vasprun = self.get_safe_vasprun(
                    file,
                    converged=False,
                    parse_dos=False,
                    parse_eigen=False,
                    parse_potcar_file=False
                )
            except FileNotFoundError:
                pass
            else:
                if vasprun is not None:
                    return vasprun.final_structure

        if try_contcar:
            try:
                file = os.path.join(run_dir, "CONTCAR")
                structure = Poscar.from_file(file).structure
                return structure
            except FileNotFoundError:
                pass

        if try_xdatcar:
            try:
                file = os.path.join(run_dir, "XDATCAR")
                structure = Xdatcar(file).structures[-1]
                return structure
            except FileNotFoundError:
                pass

        if try_poscar:
            try:
                file = os.path.join(run_dir, "POSCAR")
                structure = Poscar.from_file(file).structure
                return structure
            except FileNotFoundError:
                pass
        
        msg = (
            f"{self.get_struct_from_run_dir.__qualname__}: "
            f"Unable to extract a structure from '{run_dir}'.\n"
            "Tried files:\n"
            f"- vasprun.xml: {try_vasprun}\n"
            f"- CONTCAR: {try_contcar}\n"
            f"- XDATCAR: {try_xdatcar}\n"
            f"- POSCAR: {try_poscar}\n"
        )
        raise_or_warn(on_error, VaspParsingError, msg)
        return None


    def parse_structures(
        self,
        try_vasprun: bool = True,
        try_contcar: bool = True,
        try_xdatcar: bool = False,
        try_poscar: bool = False,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> dict[str, Structure]:
        """
        Parse all runs matching initialized indices into structures.

        Parameters
        ----------
        try_vasprun: bool
            Whether to try to get the final structure from vasprun.xml. Defaults to True.

        try_contcar: bool
            Whether to try to get the final structure from CONTCAR. Defaults to True.

        try_xdatcar: bool
            Whether to try to get the last structure from XDATCAR. Defaults to False.

        try_poscar: bool
            Whether to try to get the initial structure from POSCAR. Defaults to False.

        on_error: 'raise', 'warn', 'ignore'
            What to do when no structure could be extracted from one of the run directories.
            Defaults to "warn".

        Returns
        -------
        dict[str, Structure | None]
            Dict of run dir names and parsed Structure objects.
            Returns `None` for structures that could not be parsed.

        Raises
        ------
        `VaspParsingError` if a structure could not be extracted and `on_error`
        is set to "raise".
        """
        description = "Parsing VASP runs into structures"
        if self.workers == 0:
            struct_dict = {
                run_dir.name: self.get_struct_from_run_dir(
                    run_dir, try_vasprun, try_contcar, try_xdatcar, try_poscar, on_error
                ) for run_dir in VisualIterator(
                    self.struct_dirs,
                    desc=description,
                    unit="runs parsed",
                    percent=True
                )
            }
        else:
            get_struct_from_run_dir = ft.partial(
                self.get_struct_from_run_dir,
                try_vasprun=try_vasprun,
                try_contcar=try_contcar,
                try_xdatcar=try_xdatcar,
                try_poscar=try_poscar,
                on_error=on_error
            )
            structs = process_map(
                get_struct_from_run_dir,
                self.struct_dirs,
                max_workers=self.workers,
                chunksize=min(10, len(self.struct_dirs) // 100 + 1),
                desc=description
            )
            struct_dict = {
                run_dir.name: struct for run_dir, struct in zip(self.struct_dirs, structs)
            }

        return struct_dict


class ExtractMethod(Enum):
    FINAL_STATE = "final_state"
    CONVEX_HULL = "convex_hull"
    RELAXATION = "relaxation"
    DSOL_BANDGAP = "dsol_bandgap"


class VaspExtractor:
    """Extract data from VASP run directories."""
    def __init__(
        self,
        vasp_parser: VaspParser,
        method: str | ExtractMethod,
        summary_name: str | None = None,
        summary_key: str | None = None,
        workers: int | None = None,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> None:
        """
        Extract data from VASP run directories.

        Parameters
        ----------
        vasp_parser: VaspParser
            Initialized VaspParser object pointing to a VASP runs base directory.

        method: str | ExtractMethod
            Extraction method to use for data retrieving, depending on the kind of run to extract.

        summary_name: str, optional
            Name of the JSON file summarizing the runs results, and containing a boolean value for
            each run saying if it passed or failed the computation step filter. Avoids to load
            structures from runs that failed the filter. If not given, all run directories will be
            parsed. The file must be located inside the base diretory pointed by the VaspIO object.

        summary_key: str, optional
            Name of the key to check the filter bool value inside the summary file. If `summary_name`
            is given, this argument is mandatory.

        workers: int, optional
            Number of parallel processes to spawn for efficient file parsing. If not given,
            `tqdm.contrib.concurrent.process_map()` default is used. Pass 0 to disable
            `process_map()` and execute sequentially.

        on_error: "raise", "warn", "ignore"
            What to do in case of an error occurring while parsing VASP Δ-Sol runs data.
            Defaults to "warn".
        """
        assert isinstance(method, (str, ExtractMethod)), TypeError(
            f"'method' expected a type 'str' or 'ExtractMethod', got {type(method).__name__}."
        )
        if isinstance(method, str):
            try:
                method = ExtractMethod(method)
            except ValueError:
                supported_methods = ", ".join(sorted([method.value for method in ExtractMethod]))
                raise NotImplementedError(
                    f"Given 'method' value ({method}) is not supported. "
                    f"Supported methods are: {supported_methods}."
                )

        if summary_name is None:
            summary_path = None
        else:
            summary_path = os.path.join(vasp_parser.base_dir, summary_name)
            check_file_or_dir(summary_path, "file", allowed_formats="json")
            assert summary_key is not None, ValueError(
                "'summary_name' has been given, but not 'summary_key'. "
                "Please provide a key to match inside the summmary file."
            )

        if workers is not None:
            assert workers >= 0, ValueError(
                f"'workers' must be positive or zero, got {workers}."
            )

        self.parser = vasp_parser
        self.extractor = self._get_extractor(method, summary_path, summary_key, on_error)
        self.struct_dirs = self.parser.struct_dirs

        if workers is not None and workers == 0: # Sequential execution
            structs_data_list = list(filter(
                None, [self.extractor(struct_dir) for struct_dir in self.struct_dirs]
            ))
        else: # Multi-process execution
            structs_data_list = list(
                filter(
                    None,
                    process_map(
                        self.extractor,
                        self.struct_dirs,
                        max_workers=workers,
                        chunksize=min(10, len(self.struct_dirs) // 100 + 1),
                        desc="Extracting data from VASP output",
                    ),
                )
            )
        self._data = dict(structs_data_list)

    def get_data(self) -> dict[str, dict[str, tp.Any]]:
        """
        Get all extracted data in a dict of parsed run directories.

        Returns
        -------
        dict[str, dict[str, Any]]
            Nested dict of the form {'run_dir name': {'data_name': data}}.
            - For 'final_state' extraction method, each subdict contains
            the final structure and its total energy (in eV).

            - For 'convex_hull' extraction method, each subdict contains
            the structure ID (run_dir name), chemical composition as
            Composition object, and initial energy (in eV), with correct
            keys to pass into phase diagram construction.

            - For 'relaxation' extraction method, each subdict contains
            the initial and final structures.

            - For 'dsol_bandgap' extraction method, each subdict contains
            the unmodified structure and total energies from all sub-runs.
        """
        return self._data

    def _get_extractor(
        self,
        method: ExtractMethod,
        summary_path: PathLike | None = None,
        summary_key: str | None = None,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> Callable:
        """Initialize data extraction method."""
        match method:
            case ExtractMethod.FINAL_STATE | ExtractMethod.CONVEX_HULL | ExtractMethod.RELAXATION:
                return ft.partial(
                    self._get_run_data,
                    method=method, summary_path=summary_path, summary_key=summary_key
                )
            case ExtractMethod.DSOL_BANDGAP:
                return ft.partial(
                    self._get_dsol_bandgap_data,
                    summary_path=summary_path, summary_key=summary_key, on_error=on_error
                )
            case str():
                raise NotImplementedError(
                    f"Provided method ({method}) is not supported. "
                    "Supported methods are: 'convex_hull', 'delta_sol_init', 'delta_sol_calc'."
                )
            case _:
                raise TypeError(
                    f"'method' expected a type 'str', got {type(method).__name__}."
                )

    @staticmethod
    def _check_summary_data(
        struct_dir: PathLike, summary_file: PathLike, summary_key: str
    ) -> bool:
        """
        Check whether given structure directory is present in the summary file and
        whether it was rejected in the previous step.
        """
        with open(summary_file, "rt", encoding="utf-8") as fp:
            summary = json.load(fp)

        try:
            prev_struct_data = next(filter(
                lambda data: os.path.samefile(data["path"], struct_dir),
                summary
            ))
        except StopIteration:
            warnings.warn(
                f"Structure path '{struct_dir}' was not found in the summary file.\n"
                "You may want to pass it to the previous screening step before this one.\n"
                "This structure is assumed not viable and is ignored for this step.",
                stacklevel=2
            )
            return False

        return prev_struct_data[summary_key]

    def _get_run_data(
        self,
        struct_dir: PathLike,
        method: ExtractMethod,
        summary_path: PathLike | None = None,
        summary_key: str | None = None
    ) -> tuple[str, dict[str, tp.Any]] | None:
        """Extracts VASP data from a run directory."""
        if (
            summary_path is not None and summary_key is not None and
            not self._check_summary_data(struct_dir, summary_path, summary_key)
        ):
                return None

        vasprun = self.parser.get_safe_vasprun(
            struct_dir, parse_dos=False, parse_eigen=False, parse_potcar_file=False
        )

        if vasprun is None:
            return None

        struct_name  = os.path.basename(struct_dir)
        match method:
            case ExtractMethod.FINAL_STATE:
                struct_dict = {
                    "structure": vasprun.final_structure,
                    "final_energy": vasprun.final_energy,
                }
            case ExtractMethod.CONVEX_HULL:
                # We want data of the generated structure for convex hulls, not the relaxed one
                struct_dict = {
                    "entry_id": struct_name,
                    "composition": vasprun.initial_structure.composition,
                    "final_energy": vasprun.ionic_steps[0]["e_0_energy"]
                }
            case ExtractMethod.RELAXATION:
                struct_dict = {
                    "in_struct": vasprun.initial_structure,
                    "out_struct": vasprun.final_structure
                }
            case str():
                raise NotImplementedError(f"Method {method.value!r} not implemented.")
            case _:
                raise TypeError(
                    f"'method' expected a type 'ExtractMethod', got {type(method).__name__!r}."
                )

        return struct_name, struct_dict

    def _get_dsol_bandgap_data(
        self,
        struct_dir: PathLike,
        summary_path: PathLike | None = None,
        summary_key: str | None = None,
        on_error: tp.Literal["raise", "warn", "ignore"] = "warn"
    ) -> tuple[str, dict[str, Structure | float]] | None:
        """
        Extract the results of all Δ-Sol computations for one structure.
        """
        calc_dirs = list(filter(os.path.isdir, os.listdir(struct_dir)))
        struct_dict = {}

        for calc_dir in calc_dirs:
            if (
                summary_path is not None and summary_key is not None
                and not self._check_summary_data(struct_dir, summary_path, summary_key)
            ):
                return None

            calc_data = self._get_run_data(
                calc_dir, ExtractMethod.FINAL_STATE, summary_path, summary_key
            )
            if not calc_data:
                msg = (
                    f"{calc_dir} could not be parsed, either because the "
                    "VASP run terminated on an error, on a timeout limit "
                    "or it did not converge after the maximum ionic step was reached."
                )
                raise_or_warn(on_error, VaspParsingError, msg)
                continue

            struct_dict.setdefault("structure", calc_data[1]["structure"])
            struct_dict[calc_data[0]] = calc_data[1]["final_energy"]

        struct_name = Path(struct_dir).name

        return struct_name, struct_dict
