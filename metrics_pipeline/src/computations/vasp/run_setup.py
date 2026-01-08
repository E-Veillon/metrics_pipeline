"""Fonctions to setup and manage input data sets for VASP files."""

import typing as tp

from pymatgen.core import SiteCollection, Structure
from pymatgen.io.vasp.sets import VaspInput, VaspInputSet

from .presets import PMGRelaxSet, PMGStaticSet, DSolStaticSet
from src.utils import get_all_valence_electrons
from src.computations.local import DSolInput


U_VALUES = {
    "F": {
        "Ag": 1.5, "Co": 3.4, "Cr": 3.5, "Cu": 4.0, #"Cu": 4 -> 4.0
        "Fe": 4.0, "Mn": 3.9, "Mo": 3.5, "Nb": 1.5, #"Mo": 4.38 -> 3.5 (according to the reference)
        "Ni": 6.0, "Re": 2.0, "Ta": 2.0, "V": 3.1,  #"Ni": 6 -> 6.0, "Re": 2 -> 2.0, "Ta": 2 -> 2.0
        "W": 4.0
    },
    "O": {
        "Ag": 1.5, "Co": 3.4, "Cr": 3.5, "Cu": 4.0, #"Cu": 4 -> 4.0
        "Fe": 4.0, "Mn": 3.9, "Mo": 3.5, "Nb": 1.5, #"Mo": 4.38 -> 3.5 (according to the reference)
        "Ni": 6.0, "Re": 2.0, "Ta": 2.0, "V": 3.1,  #"Ni": 6 -> 6.0, "Re": 2 -> 2.0, "Ta": 2 -> 2.0
        "W": 4.0                          
    },
    "S": {
        "Fe": 1.9, "Mn": 2.5
    }}
"""
Values of the Hubbard U correction used in GGA + U framework, as fitted by Jain et al.

Reference:
- A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder,
Computational Materials Science, 50, 2295-2310 (2011)
"""


def _mitrelaxset_incar_corrections(n_sites: int|None = None) -> dict[str, tp.Any]:
    """
    Systematic correction for MITRelaxSet INCAR tags that do not match with 
    parameters given in the original work of Jain et al.

    Reference:
        A. Jain, G. Hautier, C.J. Moore, S.P. Ong, C.C. Fischer, T. Mueller,
        K.A. Persson, and G. Ceder, Computational Materials Science, 50, 2295-2310 (2011).

    Parameters:
        n_sites (int):  The number of sites in the structure. Used to pass explicitly
                        an EDIFF tag to pymatgen. Usually it is not necessary, as
                        the EDIFF_PER_ATOM tag is supported for the same function.
                        If not given, the correction will replace the default EDIFF
                        by an EDIFF_PER_ATOM of 5e-5 eV/atom.

    Returns:
        A dictionnary containing the INCAR tags corrections for MITRelaxSet.
    """
    corrections = {}

    if n_sites is None:
        corrections.update(
            {
                "EDIFF": None,
                "EDIFF_PER_ATOM": round(float(5e-5), 6)
            }
        )
    else:
        corrections.update({"EDIFF": round(float(5e-5) * n_sites, 6)})

    corrections.update(
        {
            "LDAUU": U_VALUES,
            "LDAUL": {
                "F": {
                    "Ag": 2, "Co": 2, "Cr": 2, "Cu": 2, "Fe": 2,
                    "Mn": 2, "Mo": 2, "Nb": 2, "Ni": 2, "Re": 2,
                    "Ta": 2, "V": 2, "W": 2
                },
                "O": {
                    "Ag": 2, "Co": 2, "Cr": 2, "Cu": 2, "Fe": 2,
                    "Mn": 2, "Mo": 2, "Nb": 2, "Ni": 2, "Re": 2,
                    "Ta": 2, "V": 2, "W": 2
                },
                "S": {
                    "Fe": 2,
                    "Mn": 2,  #"Mn": 2.5 -> 2 (quantum number l have to be an integer)
                },
            }
        }
    )
    return corrections


def _relax_set_init(
    structure: SiteCollection,
    preset: str,
    corrections: dict[str, tp.Any] | None = None,
) -> VaspInputSet:
    """Init a relaxation set of VASP input files."""
    corrections = {} if corrections is None else corrections
    incar_corrections = {}

    if preset.lower() == PMGRelaxSet.MITRELAXSET.name.lower():
        incar_corrections = _mitrelaxset_incar_corrections(structure.num_sites)

    incar_corrections.update(corrections.get("INCAR", {}))

    vasp_input_set = PMGRelaxSet.get_preset(preset)(
        structure=structure,
        user_incar_settings = incar_corrections,
        user_kpoints_settings = corrections.get("KPOINTS", {}),
        user_potcar_settings = corrections.get("POTCAR", {}),
        user_potcar_functional = corrections.get("POTCAR_FUNCTIONAL", {}),
    )

    return vasp_input_set


def _static_set_init(
    structure: SiteCollection,
    preset: str,
    nelect: float | None = None,
    corrections: dict[str, tp.Any] | None = None,
) -> VaspInputSet:
    """Init a static calculation set of VASP input files."""
    corrections = {} if corrections is None else corrections

    if preset == "DSolStaticSet":
        vasp_input_set = DSolStaticSet(
            structure=structure,
            incar_nelect=nelect,
            user_incar_settings=corrections.get("INCAR", {}),
            user_kpoints_settings=corrections.get("KPOINTS", {}),
            user_potcar_settings=corrections.get("POTCAR", {}),
            user_potcar_functional=corrections.get("POTCAR_FUNCTIONAL", {}),
        )
    else:
        vasp_input_set = PMGStaticSet.get_preset(preset)(
            structure=structure,
            user_incar_settings = corrections.get("INCAR", {}),
            user_kpoints_settings = corrections.get("KPOINTS", {}),
            user_potcar_settings = corrections.get("POTCAR", {}),
            user_potcar_functional = corrections.get("POTCAR_FUNCTIONAL", {}),
        )

    return vasp_input_set


def init_vasp_settings(
    structure: SiteCollection,
    preset: str,
    nelect: float | None = None,
    user_corrections: dict[str, tp.Any] | None = None,
) -> VaspInput:
    """
    Setup VASP inputs for a given structure using one of the pymatgen or local presets
    as base.

    Parameters:
        structure (SiteCollection):     The structure to write VASP inputs for.

        preset (str):                   The pymatgen preset to use for VASP inputs
                                        initialization. Can also be the homemade
                                        "DSolStaticSet" if Δ-Sol method by Chan et al.
                                        (2010) is used.

        nelect (float):                 Only used if DSolStaticSet is used.
                                        Sets the NELECT tag in INCAR file.
                                        Ignored if from_prev_calc is True.

        user_corrections (dict):        User defined settings. It allows to override
                                        some of the preset INCAR, KPOINTS or POTCAR
                                        settings if necessary. Defaults to None.
    """
    assert isinstance(structure, SiteCollection), TypeError(
        f"'structure' expected a type 'SiteCollection', got {type(structure).__name__}."
    )
    if preset == "DSolStaticSet":
        return _static_set_init(structure, preset, nelect, user_corrections).get_input_set()

    if PMGStaticSet.is_preset(preset):
        return _static_set_init(structure, preset, corrections=user_corrections).get_input_set()

    if PMGRelaxSet.is_preset(preset):
        return _relax_set_init(structure, preset, user_corrections).get_input_set()

    raise NotImplementedError(
        f"'preset' argument not recognized ({preset}). "
        "It must be one of the allowed pymatgen or local presets (case insensitive): "
        f"{', '.join(PMGStaticSet.names())}, {', '.join(PMGRelaxSet.names())}, "
        f"{DSolStaticSet.__name__}."
    )


def dsol_calc_init(
        structure: Structure,
        calc_index: int,
        preset: str = "DSolStaticSet",
        user_corrections: dict[str, tp.Any] | None = None,
    ) -> VaspInput:
    """
    Initializes one of the static calculations used for Δ-Sol method for one structure.

    Reference of the Δ-Sol method:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        (values in Table I)

    Parameters:
        structure (Structure):      The input structure.

        calc_index (int):           An integer corresponding to a delta-sol static calculation:
                                    0 = E(N0), 
                                    1-2 = E(N0 + n), E(N0 - n) respectively, using N*_best, 
                                    3-4 = E(N0 + n), E(N0 - n) respectively, using N*_min, 
                                    5-6 = E(N0 + n), E(N0 - n) respectively, using N*_max.

        preset (str):               A pymatgen VASP static preset, or the homemade
                                    DSolStaticSet. Defaults to DSolStaticSet.

        user_corrections (dict):    Additional corrections provided by the user in a
                                    separate .yaml file.

    Returns:
        The corresponding VaspInput object.
    """
    assert isinstance(structure, Structure), TypeError(
        f"'structure' expected a type 'Structure', got {type(structure).__name__}."
    )
    assert isinstance(calc_index, int), TypeError(
        f"'calc_index' expected a type 'int', got {type(calc_index).__name__}."
    )
    assert 0 <= calc_index <= 6, ValueError(
        f"'calc_index' must be between 0 and 6 included."
    )
    assert PMGStaticSet.is_preset(preset) or preset.lower() == "DSolStaticSet".lower(), (
        ValueError(
            f"'preset' got unsupported value {preset!r}. "
            "Supported presets (case insensitive): "
            f"{', '.join(PMGStaticSet.names() + ['DSolStaticSet'])}."
        )
    )

    nb_val_elec = get_all_valence_electrons(structure)
    run_set = init_vasp_settings(structure, preset, user_corrections=user_corrections)

    # Search for the right N* parameter to use with respect to the functional
    pot_func = run_set.get("POTCAR_FUNCTIONAL", "PBE")

    n_ratio = DSolInput(
        structure=structure,
        calc_type=calc_index,
        functional=pot_func
    ).n_ratio

    nelect = nb_val_elec + n_ratio if calc_index % 2 == 1 else nb_val_elec - n_ratio

    if preset == "DSolStaticSet":
        run_set = init_vasp_settings(
            structure, preset, nelect=nelect, user_corrections=user_corrections
        )

    else:
        run_dict = run_set.as_dict()
        run_dict["INCAR"].update({"NELECT": nelect})
        run_set = VaspInput.from_dict(run_dict)

    return run_set
