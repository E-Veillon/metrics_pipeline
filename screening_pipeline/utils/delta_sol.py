"""Functions that are specific to Δ-Sol method application."""


import os
from typing import Tuple, Dict, Union, Literal, Any, Optional
from dataclasses import dataclass
from tqdm.contrib.concurrent import process_map

# PYTHON MATERIALS GENOMICS
from pymatgen.core import SiteCollection, Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import VaspInput
from pymatgen.io.vasp.sets import MPRelaxSet

# LOCAL IMPORTS
from periodic_table import get_all_valence_electrons
from file_io import yaml_loader
from vasp_io import vasp_static_settings
from custom_types import PathLike, PMGStaticSetType, PMGStaticSet
from fitted_values import EL_PER_XC_VOL


########################################


@dataclass
class DeltaSolStaticSet(MPRelaxSet):
    """
    Initialize VASP input files for Δ-Sol method computations using
    PBE_54_W_HASH pymatgen set of POTCAR files. Parameters are as 
    described in Δ-Sol method original work by Chan et al. in 2010.
    DFT+U corrections are used as proposed by Jain et al. in 2011.

    References:
        - M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
        (Ref 32 in screening_pipeline/Bibliography)

        - A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller, 
        K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).
        (Ref 14 in screening_pipeline/Bibliography)

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
    base_path = os.path.dirname(os.path.dirname(__file__))
    path = os.path.join(base_path, "config", "DeltaSolStaticSet.yaml")
    CONFIG = yaml_loader(path, on_error="raise")

    def __init__(
            self,
            structure: Structure|None = None,
            incar_nelect: float|None = None,
            **kwargs
        ) -> None:
        """DeltaSolStaticSet init."""
        super().__init__(structure, **kwargs)

        if incar_nelect is None:
            try:
                incar_nelect = get_all_valence_electrons(structure)
            except TypeError as exc:
                raise ValueError("Either structure or incar_nelect must be given.") from exc

        self.incar_nelect = incar_nelect

    @property
    def incar_updates(self) -> Dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: Dict[str, Any] = {"MAGMOM": None, "NELECT": self.incar_nelect}
        return updates


########################################


def get_dsol_struct_dir(
        path: PathLike, task_id: int, with_uncertainties: bool = False
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
    except StopIteration as exc:
        raise ValueError(
            "No structure directory found with index corresponding to given 'task-id' argument.\n"
            f"Searched directory: {path}\n"
            f"Given task-id argument: {task_id}\n"
            f"Corresponding structure index: {struct_idx}\n"
            f"Corresponding calculation ID: {calc_idx}\n"
            "(0 = E(N0), 1-2 = E(N0 +/- n(best)), 3-4 = E(N0 +/- n(min)), 5-6 = E(N0 +/- n(max)))."
        ) from exc

    struct_path = os.path.join(path, struct_dir)

    return struct_path, struct_idx, calc_idx


########################################


def calc_idx_to_dir_name(struct_dir_name: str, calc_index: int) -> str:
    """Maps calculation index to corresponding calculation name."""
    match calc_index:
        case 0:
            calc_type = "neutral"
        case 1:
            calc_type = "best_plus"
        case 2:
            calc_type = "best_minus"
        case 3:
            calc_type = "min_plus"
        case 4:
            calc_type = "min_minus"
        case 5:
            calc_type = "max_plus"
        case 6:
            calc_type = "max_minus"
        case int():
            raise ValueError(
                f"Only int from 0 to 6 supported, got {calc_index}"
            )
        case _:
            raise TypeError(
                "'calc_index' expected a type 'int', "
                f"got '{type(calc_index)}' instead."
            )
    return "_".join((struct_dir_name, calc_type))


########################################


def get_dsol_n_ratio(
        structure: SiteCollection,
        dft_functional: Literal['LDA','PBE','AM05'] = 'PBE',
        n_star_type: Literal['MIN', 'BEST', 'MAX'] = 'BEST'
    ) -> float:
    '''
    Computes n = N0/N* the electron ratio to add or remove from 
    the structure in the Δ-Sol method developped by Chan et al.

    Reference:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (reference 32 in screening_pipeline/Bibliography)
    '''
    val_elec_type = 'sp'

    for elt in structure.elements:
        if elt.block in {'s', 'p'}:
            continue
        if elt.block == 'd':
            val_elec_type = 'spd'
            break
        if elt.block == 'f':
            raise NotImplementedError(
                'f-block elements are not supported in Δ-Sol method.'
            )
        raise ValueError(
            'Something is wrong with this function or Element objects "block" property.'
        )

    n_0        = get_all_valence_electrons(structure)
    value_name = '_'.join((dft_functional, val_elec_type))
    n_star     = EL_PER_XC_VOL[n_star_type][value_name]
    n          = float(n_0) / float(n_star)

    return n


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


def _match_calc_index(calc_index: int) -> str|None:
    match calc_index:
        case 0:
            return None
        case 1|2:
            return "BEST"
        case 3|4:
            return "MIN"
        case 5|6:
            return "MAX"
        case int():
            raise ValueError("calc_index must be between 0 and 6 included.")
        case _:
            raise TypeError(f"Expected 'int' type, got '{type(calc_index)}' type instead.")


########################################


def dsol_calc_init(
        structure: Structure,
        calc_index: int,
        preset: PMGStaticSetType|"DeltaSolStaticSet" = "DeltaSolStaticSet",
        user_corrections: Optional[Dict[str, Any]] = None,
    ) -> VaspInput:
    """
    Initializes one of the static calculations used for Δ-Sol method for one structure.

    Reference of the Δ-Sol method:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (reference 32 in screening_pipeline/Bibliography, values in Table I)

    Parameters:
        structure (Structure):      The input structure.

        calc_index (int):           An integer corresponding to a delta-sol static calculation:
                                    0 = E(N0), 
                                    1-2 = E(N0 + n), E(N0 - n) respectively, using N*_best, 
                                    3-4 = E(N0 + n), E(N0 - n) respectively, using N*_min, 
                                    5-6 = E(N0 + n), E(N0 - n) respectively, using N*_max.

        preset (str):               A pymatgen VASP static preset, or the homemade
                                    DeltaSolStaticSet. Defaults to DeltaSolStaticSet.

        user_corrections (dict):    Additional corrections provided by the user in a
                                    separate .yaml file.

    Returns:
        The corresponding VaspInput object.
    """
    if not isinstance(structure, Structure):
        raise TypeError(
            "'structure' argument expected a type 'pymatgen.core.structure.Structure', "
            f"got '{type(structure)}' instead."
        )
    if not isinstance(calc_index, int):
        raise TypeError(
            "'calc_index' argument expected a type 'int', "
            f"got '{type(calc_index)}' instead."
        )
    if not 0 <= calc_index <= 6:
        raise ValueError(
            "'calc_index' argument value must be between 0 and 6 included, "
            f"got {calc_index}."
        )
    if not preset in PMGStaticSet or preset == "DeltaSolStaticSet":
        raise ValueError(
            "'preset' argument value is not a supported preset. "
            "Supported presets are:\n"
            f"{PMGStaticSet + set(('DeltaSolStaticSet',))}\n"
            f"'preset' got value '{preset}' instead."
        )
    if not isinstance(user_corrections, Dict) and user_corrections is not None:
        raise TypeError(
            "'user_corrections' expected a type 'dict', "
            f"got '{type(user_corrections)}' instead."
        )

    nb_val_elec = get_all_valence_electrons(structure)
    run_set = vasp_static_settings(structure, preset, user_corrections=user_corrections)

    # Search for the right N* parameter to use with respect to the functional
    pot_func = run_set.get("POTCAR_FUNCTIONAL", "PBE")
    dsol_functional = _match_dft_functional(pot_func)
    n_star_type = _match_calc_index(calc_index)

    if n_star_type is not None:
        n_ratio = get_dsol_n_ratio(
            structure=structure,
            dft_functional=dsol_functional,
            n_star_type=n_star_type
        )
        nelect = nb_val_elec + n_ratio if calc_index % 2 == 1 else nb_val_elec - n_ratio

        if preset == "DeltaSolStaticSet":
            run_set = vasp_static_settings(
                structure, preset, nelect=nelect, user_corrections=user_corrections
            )
    else:
        nelect = nb_val_elec

    if preset != "DeltaSolStaticSet":
        run_dict = run_set.as_dict()
        run_dict["INCAR"].update({"NELECT": nelect})
        run_set = VaspInput.from_dict(run_dict)

    return run_set


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

    def has_str_key(dct: Dict, key: str) -> bool:
        return dct.get(key) is not None

    data_keys = ("name","functional","structure","E_N0","E_N0_plus_n_best","E_N0_minus_n_best")
    supp_keys = ("E_N0_plus_n_min","E_N0_minus_n_min","E_N0_plus_n_max","E_N0_minus_n_max")

    assert isinstance(data, dict)
    assert all(has_str_key(data, key) for key in data_keys)

    n_ratio_best = get_dsol_n_ratio(
        data["structure"], dft_functional=data["functional"], n_star_type="BEST"
    )

    # E_FG = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n -> Δ-Sol band gap
    # (Ref 32 in screening_pipeline/Bibliography))
    e_diff_best = data["E_N0_plus_n_best"] + data["E_N0_minus_n_best"] - 2*data["E_N0"]
    e_bg_best = e_diff_best / n_ratio_best

    if not all(has_str_key(data, key) for key in supp_keys):
        # If data do not have uncertainty keys, return here
        return (data["name"], e_bg_best)

    n_ratio_min = get_dsol_n_ratio(
        data["structure"], dft_functional=data["functional"], n_star_type="MIN"
    )
    n_ratio_max = get_dsol_n_ratio(
        data["structure"], dft_functional=data["functional"], n_star_type="MAX"
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
        workers: int = 1
    ) -> Dict[str, float]:
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

        workers (int):              The number of parallel processes to spawn. Defaults to 1.

    Returns: Dict[str, float]:
        Dict of Band gap values associated with the original structure directory name.
    """

    assert isinstance(bg_data, dict)
    assert dft_functional in {"LDA", "PBE", "AM05"}
    assert isinstance(with_uncertainties, bool)
    assert isinstance(workers, int) and workers >= 1

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
                'E_band_gap': tup[1]
            } for tup in e_band_gaps
        }

    else:
        e_band_gaps = {
            tup[0]: {
                'E_band_gap': tup[1],
                'E_band_gap_min': tup[2],
                'E_band_gap_max': tup[3]
            } for tup in e_band_gaps
        }

    return e_band_gaps


########################################


if __name__ == "__main__":
    # Test for the DeltaSolStaticSet class. You may have to change the given path
    # to one pointing at a valid CIF structure file for it to work properly.
    # You'll also need to set PMG_VASP_PSP_DIR for POTCAR files in .pmgrc.yaml.
    test_path = "/home/elohan/screening-pipeline/screening_pipeline/_benchmarks/TiO2.cif"
    with open(test_path, "rt", encoding="utf-8") as test_file:
        struct = CifParser(test_file).parse_structures()[0]
    dset = DeltaSolStaticSet(struct).get_input_set()
    with open(
        os.path.join(DeltaSolStaticSet.base_path, "DeltaVaspInput.txt"),
        mode="wt", encoding="utf-8"
    ) as out:
        out.write(str(dset))
    print(dset.incar_nelect)
    dset.incar_nelect = 58
    print(dset.incar_nelect)

    # Unit test for _match_calc_index().
    print("Wanted output:")
    print("None\nBEST BEST\nMIN MIN\nMAX MAX")
    print("Actual output:")
    print(_match_calc_index(0))
    print(_match_calc_index(1), _match_calc_index(2))
    print(_match_calc_index(3), _match_calc_index(4))
    print(_match_calc_index(5), _match_calc_index(6))
