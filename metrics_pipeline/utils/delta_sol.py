#!/usr/bin/python
"""Functions that are specific to Δ-Sol method application."""


import os
from typing import Tuple, Dict, Union, Literal, Any
from dataclasses import dataclass
from enum import Enum
from tqdm.contrib.concurrent import process_map


# PYTHON MATERIALS GENOMICS
from pymatgen.core import SiteCollection, Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.sets import MPRelaxSet

# LOCAL IMPORTS
try:
    from .common_asserts import check_type, check_num_value
    from .custom_types import PathLike
    from .file_io import CONFIGPATH, check_file_or_dir, yaml_loader
    from .fitted_values import EL_PER_XC_VOL
    from .periodic_table import get_all_valence_electrons
except ImportError as exc:
    from common_asserts import check_type, check_num_value
    from custom_types import PathLike
    from file_io import CONFIGPATH, check_file_or_dir, yaml_loader
    from fitted_values import EL_PER_XC_VOL
    from periodic_table import get_all_valence_electrons


########################################


@dataclass
class DSolStaticSet(MPRelaxSet):
    """
    Initialize VASP input files for Δ-Sol method computations using
    PBE_54_W_HASH pymatgen set of POTCAR files. Parameters are as 
    described in Δ-Sol method original work by Chan et al. in 2010.
    DFT+U corrections are used as proposed by Jain et al. in 2011.

    References:
        - M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).

        - A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller, 
        K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).

    Args:
        structure (Structure):  The Structure to create inputs for. If None, the input
                                set is initialized without a Structure but one must be
                                set separately before the inputs are generated.

        incar_nelect (float):   The number of electrons to put in the NELECT INCAR tag.
                                In Δ-Sol, several computations with distinct number of 
                                electrons are done, this is a convenient arg to set that.
                                If not given, infers the Δ-Sol N0 electrons calculation 
                                from the given structure.

        **kwargs:               kwargs supported by DictSet.
    
    Raises: ValueError if neither structure nor nelect are given at instanciation time.
    """
    CONFIG = yaml_loader(os.path.join(CONFIGPATH, "DSolStaticSet.yaml"), on_error="raise")

    def __init__(
            self,
            structure: Structure|None = None,
            incar_nelect: float|None = None,
            **kwargs
        ) -> None:
        """DSolStaticSet init."""
        super().__init__(structure, **kwargs) # type: ignore

        if incar_nelect is None:
            try:
                incar_nelect = get_all_valence_electrons(structure)
            except TypeError as e:
                raise ValueError("Either structure or incar_nelect must be given.") from e

        self.incar_nelect = incar_nelect

    @property
    def incar_updates(self) -> Dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: Dict[str, Any] = {"MAGMOM": None, "NELECT": self.incar_nelect}
        return updates


########################################


class DSolCalc(Enum):
    """Enum class to store possible Delta-Sol calculations."""
    NEUTRAL = 0
    BEST_PLUS = 1
    BEST_MINUS = 2
    MIN_PLUS = 3
    MIN_MINUS = 4
    MAX_PLUS = 5
    MAX_MINUS = 6


########################################


def get_dsol_struct_dir(
    path: PathLike, task_id: int, with_uncertainties: bool = False #type: ignore
) -> Tuple[str, int, int]:
    """
    Determines structure index and calculation ID from the task ID,
    then finds corresponding structure directory.

    Parameters:
        path (Path|str):            Where to search for the structure directory.

        task_id (int):              The ID of the task in the job array.

        with_uncertainties (bool):  Whether to take into account uncertainty calculations.
                                    Defaults to False.

    Returns:
        Tuple[str, int, int]:
        path to the structure directory, structure index and calculation ID.
    """
    check_file_or_dir(path, "dir")
    check_type(task_id, "task_id", (int,))
    check_num_value(task_id, "task_id", ">=", 0)
    check_type(with_uncertainties, "with_uncertainties", (bool,))

    tasks_per_struct = 7 if with_uncertainties else 3
    struct_idx = task_id // tasks_per_struct
    calc_idx = task_id % tasks_per_struct

    try:
        struct_dir = next(
            filter(
            lambda dirname: os.path.isdir(dirname) and dirname.startswith(f"{struct_idx}_"),
            os.listdir(path)
            )
        )
    except StopIteration as e:
        raise ValueError(
            "No structure directory found with index corresponding to given 'task-id' argument.\n"
            f"Searched directory: {path}\n"
            f"Given task-id argument: {task_id}\n"
            f"Corresponding structure index: {struct_idx}\n"
            f"Corresponding calculation ID: {calc_idx}\n"
            "(0 = E(N0), 1-2 = E(N0 +/- n(best)), 3-4 = E(N0 +/- n(min)), 5-6 = E(N0 +/- n(max)))."
        ) from e

    struct_path = os.path.join(path, struct_dir)

    return struct_path, struct_idx, calc_idx


########################################


def _get_dsol_calc(calc_idx: int) -> str:
    for calc in DSolCalc:
        if calc_idx == calc.value:
            return calc.name
    check_num_value(calc_idx, "calc_idx", ">=", 0)
    check_num_value(calc_idx, "calc_idx", "<=", 6)
    raise NotImplementedError


########################################


def calc_idx_to_dir_name(struct_dir_name: str, calc_index: int) -> str:
    """Maps calculation index to corresponding calculation name."""
    check_type(struct_dir_name, "struct_dir_name", (str,))
    check_type(calc_index, "calc_index", (int,))

    calc_type = _get_dsol_calc(calc_index)
    return "_".join((struct_dir_name, calc_type.lower()))


########################################


def _match_n_star_idx(n_star_idx: int) -> str:
    """Get delta-sol N* type name (e.g. "BEST")."""
    n_star_type = _get_dsol_calc(n_star_idx)
    return n_star_type.split(sep="_")[0]


########################################


def _match_dft_functional(
    functional: str
) -> Union[Literal["LDA"], Literal["PBE"], Literal["AM05"]]:
    """
    Verify if given DFT functional is compatible with Δ-Sol method,
    and returns corresponding N* DFT category if it is the case,
    i.e. "LDA", "PBE", or "AM05".
    """
    if "LDA" in functional:
        return "LDA"
    if "PBE" in functional:
        return "PBE"
    if "AM05" in functional:
        return "AM05"
    raise NotImplementedError(
        "Provided POTCAR functional is not implemented for Δ-Sol method.\n"
        "Recognized functionals are 'LDA', 'PBE', and 'AM05'."
    )


########################################


def _match_orbital_types(structure: SiteCollection) -> Union[Literal["sp"], Literal["spd"]]:
    """Get valence orbital types in the structure for delta-Sol."""
    for elt in structure.elements:
        match elt.block:
            case "s"|"p":
                continue
            case "d":
                return "spd"
            case "f":
                raise NotImplementedError(
                    "f-block elements are not supported in Δ-Sol method."
                )
    return "sp"


########################################


def get_dsol_n_ratio(
        structure: SiteCollection,
        dft_functional: str = "PBE",
        n_star_idx: int = 0
    ) -> float:
    """
    Computes n = N0/N* the electron ratio to add or remove from 
    the structure in the Δ-Sol method developped by Chan et al.

    Reference:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
    """
    check_type(structure, "structure", (SiteCollection,))
    check_type(n_star_idx, "n_star_idx", (int,))
    if not 0 <= n_star_idx <= 6:
        check_num_value(n_star_idx, "n_star_idx", ">=", 0)
        check_num_value(n_star_idx, "n_star_idx", "<=", 6)

    # n_star_idx = 0 => n_star_type = "NEUTRAL" => n = 0.0 electron
    if n_star_idx == 0:
        return 0.0

    value_name = "_".join(
        (
            _match_dft_functional(dft_functional),
            _match_orbital_types(structure)
        )
    )
    n_0    = get_all_valence_electrons(structure)
    n_star = EL_PER_XC_VOL[_match_n_star_idx(n_star_idx)][value_name]
    n      = float(n_0) / float(n_star)

    return n


########################################


def get_dsol_band_gap(data: dict) -> Union[Tuple[str, float], Tuple[str, float, float, float]]:
    """
    Calculate Δ-Sol band gap value of a structure, provided a dict containing all necessary data.

    Parameters:
        data (dict): A dict containing at least following data about a structure:
                        - Its name,
                        - The DFT functional used for the calculation, 
                        - The Structure object, 
                        - Its total energy with N0 electrons, 
                        - Its total energy with N0 + n(best) electrons, 
                        - Its total energy with N0 - n(best) electrons.

                     It can also contain data for uncertainty calculations:
                        - Total energy with N0 + n(min) electrons,
                        - Total energy with N0 - n(min) electrons,
                        - Total energy with N0 + n(max) electrons,
                        - Total energy with N0 - n(max) electrons.

    Returns:
        The name of the structure and its band gap value(s).
    """
    check_type(data, "data", (Dict,))

    def has_str_key(dct: Dict, key: str) -> bool:
        return dct.get(key) is not None

    data_keys = ("name","functional","structure","E_N0","E_N0_plus_n_best","E_N0_minus_n_best")
    supp_keys = ("E_N0_plus_n_min","E_N0_minus_n_min","E_N0_plus_n_max","E_N0_minus_n_max")

    assert all(has_str_key(data, key) for key in data_keys)

    n_ratio_best = get_dsol_n_ratio(
        data["structure"], dft_functional=data["functional"], n_star_idx=1
    )

    # E_FG = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n -> Δ-Sol band gap
    # (Reference of the method: M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010))
    e_diff_best = data["E_N0_plus_n_best"] + data["E_N0_minus_n_best"] - 2*data["E_N0"]
    e_bg_best = e_diff_best / n_ratio_best

    if not all(has_str_key(data, key) for key in supp_keys):
        # If data do not have uncertainty keys, return here
        return (data["name"], e_bg_best)

    n_ratio_min = get_dsol_n_ratio(
        data["structure"], dft_functional=data["functional"], n_star_idx=3
    )
    n_ratio_max = get_dsol_n_ratio(
        data["structure"], dft_functional=data["functional"], n_star_idx=5
    )
    e_diff_min = data["E_N0_plus_n_min"] + data["E_N0_minus_n_min"] - 2*data["E_N0"]
    e_bg_min = e_diff_min / n_ratio_min

    e_diff_max = data["E_N0_plus_n_max"] + data["E_N0_minus_n_max"] - 2*data["E_N0"]
    e_bg_max = e_diff_max / n_ratio_max

    return (data["name"], e_bg_best, e_bg_min, e_bg_max)


########################################


def batch_get_dsol_band_gaps(
        bg_data: dict,
        dft_functional: Literal["LDA","PBE","AM05"] = "PBE",
        with_uncertainties: bool = False,
        workers: int | None = None
    ) -> dict[str, dict[str, float]]:
    """
    Calculate Δ-Sol band gap value for every structure in a batch
    from their data, as provided by extract_vasp_data_for_delta_sol
    function applied on the 3 energy calculations.

    Parameters:
        bg_data (dict):             Dict containing structures data extracted
                                    from previous VASP static calculations.

        dft_functional (str):       The type of functional used for static calculations.
                                    Supported functionals are "LDA", "PBE", and "AM05".
                                    Defaults to "PBE".

        with_uncertainties (bool):  Whether to include uncertainty calculations data
                                    in the results. Defaults to False.

        workers (int):              Number of parallel processes to spawn.
                                    If not given, tqdm.contrib.concurrent.process_map
                                    default is used.

    Returns:
        Dict of Band gap values associated with the original structure directory name.
    """
    check_type(bg_data, "bg_data", (Dict,))
    assert dft_functional in {"LDA", "PBE", "AM05"}
    check_type(with_uncertainties, "with_uncertainties", (bool,))
    if workers is not None:
        check_type(workers, "workers", (int,))
        check_num_value(workers, "workers", ">", 0)

    if not with_uncertainties:
        final_energies = {
            name: {
                "name": name,
                "functional": dft_functional,
                "structure": data["structure"],
                "E_N0": data[name + "_neutral"],
                "E_N0_plus_n_best": data[name + "_best_plus"],
                "E_N0_minus_n_best": data[name + "_best_minus"],
            } for name, data in bg_data.items()
        }

    else:
        final_energies = {
            name: {
                "name": name,
                "functional": dft_functional,
                "structure": data["structure"],
                "E_N0": data[name + "_neutral"],
                "E_N0_plus_n_best": data[name + "_best_plus"],
                "E_N0_minus_n_best": data[name + "_best_minus"],
                "E_N0_plus_n_min": data[name + "_min_plus"],
                "E_N0_minus_n_min": data[name + "_min_minus"],
                "E_N0_plus_n_max": data[name + "_max_plus"],
                "E_N0_minus_n_max": data[name + "_max_minus"],
            } for name, data in bg_data.items()
        }

    nbr_structs = len(final_energies)
    chunksize = min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1
    data_list = list(final_energies.values())

    e_band_gaps = list(process_map(
        get_dsol_band_gap,
        data_list,
        max_workers=workers,
        chunksize=chunksize
    ))

    if not with_uncertainties:
        e_band_gaps = {
            tup[0]: {
                "E_band_gap": tup[1]
            } for tup in e_band_gaps
        }

    else:
        e_band_gaps = {
            tup[0]: {
                "E_band_gap": tup[1],
                "E_band_gap_min": tup[2],
                "E_band_gap_max": tup[3]
            } for tup in e_band_gaps
        }

    return e_band_gaps


########################################


if __name__ == "__main__":
    # Test for the DSolStaticSet class. You may have to change the given path
    # to one pointing at a valid CIF structure file for it to work properly.
    # You'll also need to set PMG_VASP_PSP_DIR for POTCAR files in .pmgrc.yaml.
    file_path = ""
    PATHTEST = os.path.join(os.path.expanduser("~"), file_path)
    with open(PATHTEST, "rt", encoding="utf-8") as test_file:
        struct = CifParser(test_file).parse_structures()[0] # type: ignore
    dset = DSolStaticSet(struct).get_input_set()
    with open(
        os.path.join(os.path.dirname(PATHTEST), "DeltaVaspInput.txt"),
        mode="wt", encoding="utf-8"
    ) as out:
        out.write(str(dset))
    print(dset.incar_nelect) # type: ignore
    dset.incar_nelect = 58 # type: ignore
    print(dset.incar_nelect) # type: ignore

    # Unit test for _match_calc_index().
    print("Wanted output:")
    print("None\nBEST BEST\nMIN MIN\nMAX MAX")
    print("Actual output:")
    print(_match_n_star_idx(0))
    print(_match_n_star_idx(1), _match_n_star_idx(2))
    print(_match_n_star_idx(3), _match_n_star_idx(4))
    print(_match_n_star_idx(5), _match_n_star_idx(6))
