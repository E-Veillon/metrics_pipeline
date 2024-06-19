"""Functions that are specific to Δ-Sol method application."""

from typing import Tuple, Dict, Union, Literal

from tqdm.contrib.concurrent import process_map

from pymatgen.core import SiteCollection

from screening_pipeline.utils import get_all_valence_electrons


########################################


def get_delta_sol_n_ratio(
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
    from screening_pipeline.utils import EL_PER_XC_VOL
    
    val_elec_type = 'sp'

    for elt in structure.elements:
        if elt.block == 's' or elt.block == 'p':
            continue
        elif elt.block == 'd': 
            val_elec_type = 'spd'
            break
        elif elt.block == 'f':
            raise NotImplementedError(
                'f-block elements are not supported in Δ-Sol method.'
            )
        else:
            raise ValueError(
                'Something is wrong with this function or Element objects "block" property.'
            )
    
    N_0        = get_all_valence_electrons(structure)
    value_name = '_'.join((dft_functional, val_elec_type))
    N_star     = EL_PER_XC_VOL[n_star_type][value_name]
    n          = float(N_0) / float(N_star)

    return n


########################################


def calculate_delta_sol_band_gap(data: dict) -> Union[Tuple[str, float], Tuple[str, float, float, float]]:
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
    assert all([has_str_key(data, key) for key in data_keys])

    n_ratio_best = get_delta_sol_n_ratio(data["structure"], dft_functional=data["functional"], n_star_type="BEST")

    # E_FG = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n -> Δ-Sol band gap 
    # (Ref 32 in screening_pipeline/Bibliography))
    E_band_gap = (data["E_N0_plus_n_best"] + data["E_N0_minus_n_best"] - 2*data["E_N0"])/n_ratio_best

    if not all([has_str_key(data, key) for key in supp_keys]): # data do not have uncertainty keys
        return (data["name"], E_band_gap)

    n_ratio_min = get_delta_sol_n_ratio(data["structure"], dft_functional=data["functional"], n_star_type="MIN")
    n_ratio_max = get_delta_sol_n_ratio(data["structure"], dft_functional=data["functional"], n_star_type="MAX")

    E_band_gap_min = (data["E_N0_plus_n_min"] + data["E_N0_minus_n_min"] - 2*data["E_N0"])/n_ratio_min
    E_band_gap_max = (data["E_N0_plus_n_max"] + data["E_N0_minus_n_max"] - 2*data["E_N0"])/n_ratio_max

    return (data["name"], E_band_gap, E_band_gap_min, E_band_gap_max)


########################################


def batch_calculate_delta_sol_band_gaps(
        bg_data: dict, 
        dft_functional: Literal["LDA","PBE","AM05"] = "PBE", 
        with_uncertainties: bool = False, 
        workers: int = 1, 
        /
    ) -> Dict[str, float]:
    """
    Calculate Δ-Sol band gap value for every structure in a batch from their data,
    as provided by extract_vasp_data_for_delta_sol function applied on the 3 energy calculations.

    Parameters:
        bg_data (dict):               Dict containing structures data extracted from previous VASP static calculations.

        dft_functional (str):         The type of functional used for static calculations. Supported functionals are
                                      "LDA", "PBE", and "AM05". Defaults to "PBE".

        with_uncertainties (bool):    Whether to include uncertainty calculations data in the results.
                                      Defaults to False.

        workers (int):                The number of parallel processes to spawn. Defaults to 1.

    Returns:
        Dict[str, float]: Dict of Band gap values associated with the original structure directory name.
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

    E_band_gaps = list(process_map(
        calculate_delta_sol_band_gap, 
        data_list, 
        max_workers=workers, 
        chunksize=chunksize
    ))

    if not with_uncertainties:
        E_band_gaps = {
            tup[0]: {
                'E_band_gap': tup[1]
            } for tup in E_band_gaps
        }

    else:
        E_band_gaps = {
            tup[0]: {
                'E_band_gap': tup[1], 
                'E_band_gap_min': tup[2], 
                'E_band_gap_max': tup[3]
            } for tup in E_band_gaps
        }

    return E_band_gaps


########################################