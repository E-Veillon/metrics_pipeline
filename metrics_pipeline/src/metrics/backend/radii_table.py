"""Reference atomic radii tables for Viability metric computation in picometers."""

# Label for values reported in the Wikipedia article but not bound to a source reference
NOREF = 1

##### EXPERIMENTAL MEASUREMENTS #####

# Inferior limit for unknown values of radius
UNKNOWN = 25

# Standard deviation reported by J. C. Slater in his "Atomic Radii in Crystals" article published in 1964.
slater_std = 12

slater_radii_table_pm_1 = {
    # 1s
    "H": 25, "He": 120*NOREF,
    # 2s
    "Li": 145, "Be": 105,
    # 2p
    "B": 85, "C": 70, "N": 65, "O": 60, "F": 50, "Ne": 60*NOREF,
    # 3s
    "Na": 180, "Mg": 150,
    # 3p
    "Al": 125, "Si": 110, "P": 100, "S": 100, "Cl": 100, "Ar": 71*NOREF,
    # 4s
    "K": 220, "Ca": 180,
    # 3d
    "Sc": 160, "Ti": 140, "V": 135, "Cr": 140, "Mn": 140,
    "Fe": 140, "Co": 135, "Ni": 135, "Cu": 135, "Zn": 135,
    # 4p
    "Ga": 130, "Ge": 125, "As": 115, "Se": 115, "Br": 115, "Kr": UNKNOWN,
    # 5s
    "Rb": 235, "Sr": 200,
    # 4d
    "Y": 180, "Zr": 155, "Nb": 145, "Mo": 145, "Tc": 135,
    "Ru": 130, "Rh": 135, "Pd": 140, "Ag": 160, "Cd": 155,
    # 5p
    "In": 155, "Sn": 145, "Sb": 145, "Te": 140, "I": 140, "Xe": UNKNOWN,
    # 6s
    "Cs": 260, "Ba": 215,
    # 4f
    "La": 195, "Ce": 185, "Pr": 185, "Nd": 185, "Pm": 185, "Sm": 185, "Eu": 185,
    "Gd": 180, "Tb": 175, "Dy": 175, "Ho": 175, "Er": 175, "Tm": 175, "Yb": 175,
    # 5d
    "Lu": 175, "Hf": 155, "Ta": 145, "W": 135, "Re": 135,
    "Os": 130, "Ir": 135, "Pt": 135, "Au": 135, "Hg": 150,
    # 6p
    "Tl": 190, "Pb": 180*NOREF, "Bi": 160, "Po": 190, "At": UNKNOWN, "Rn": UNKNOWN,
    # 7s
    "Fr": UNKNOWN, "Ra": 215,
    # 5f
    "Ac": 195, "Th": 180, "Pa": 180, "U": 175, "Np": 175, "Pu": 175, "Am": 175,
    "Cm": 176*NOREF, "Bk": UNKNOWN, "Cf": UNKNOWN, "Es": UNKNOWN, "Fm": UNKNOWN, "Md": UNKNOWN, "No": UNKNOWN,
    # 6d
    "Lr": UNKNOWN, "Rf": UNKNOWN, "Db": UNKNOWN, "Sg": UNKNOWN, "Bh": UNKNOWN,
    "Hs": UNKNOWN, "Mt": UNKNOWN, "Ds": UNKNOWN, "Rg": UNKNOWN, "Cn": UNKNOWN,
    # 7p
    "Nh": UNKNOWN, "Fl": UNKNOWN, "Mc": UNKNOWN, "Lv": UNKNOWN, "Ts": UNKNOWN, "Og": UNKNOWN
}
"""
Wikipedia: "Rayons atomiques des éléments (page de données)" - Table 1, column "valeurs empiriques" (20/02/2025).

Reference :
J. C. Slater, « Atomic Radii in Crystals »,
The Journal of Chemical Physics, vol. 41, no 10, 1964, p. 3199–3204
"""

slater_radii_table_pm_2 = {
    # 1s
    "H": 25, "He": UNKNOWN,
    # 2s
    "Li": 145, "Be": 105,
    # 2p
    "B": 95, "C": 85, "N": 85, "O": 90, "F": 50, "Ne": UNKNOWN,
    # 3s
    "Na": 180, "Mg": 150,
    # 3p
    "Al": 125, "Si": 110, "P": 100, "S": 100, "Cl": 100, "Ar": UNKNOWN,
    # 4s
    "K": 220, "Ca": 180,
    # 3d
    "Sc": 160, "Ti": 140, "V": 135, "Cr": 140, "Mn": 140,
    "Fe": 140, "Co": 135, "Ni": 135, "Cu": 135, "Zn": 135,
    # 4p
    "Ga": 130, "Ge": 125, "As": 115, "Se": 115, "Br": 115, "Kr": UNKNOWN,
    # 5s
    "Rb": 265, "Sr": 200,
    # 4d
    "Y": 180, "Zr": 155, "Nb": 145, "Mo": 145, "Tc": 135,
    "Ru": 130, "Rh": 135, "Pd": 140, "Ag": 160, "Cd": 155,
    # 5p
    "In": 155, "Sn": 145, "Sb": 145, "Te": 140, "I": 140, "Xe": UNKNOWN,
    # 6s
    "Cs": 260, "Ba": 215,
    # 4f
    "La": 195, "Ce": 185, "Pr": 185, "Nd": 185, "Pm": 185, "Sm": 185, "Eu": 185,
    "Gd": 180, "Tb": 175, "Dy": 175, "Ho": 175, "Er": 175, "Tm": 175, "Yb": 175,
    # 5d
    "Lu": 175, "Hf": 155, "Ta": 145, "W": 135, "Re": 135,
    "Os": 130, "Ir": 135, "Pt": 135, "Au": 135, "Hg": 150,
    # 6p
    "Tl": 190, "Pb": 180, "Bi": 160, "Po": 190, "At": 180, "Rn": UNKNOWN,
    # 7s
    "Fr": 280, "Ra": 285,
    # 5f
    "Ac": 195, "Th": 180, "Pa": 180, "U": 175, "Np": 175, "Pu": 175, "Am": 175,
    "Cm": UNKNOWN, "Bk": UNKNOWN, "Cf": UNKNOWN, "Es": UNKNOWN, "Fm": UNKNOWN, "Md": UNKNOWN, "No": UNKNOWN,
    # 6d
    "Lr": UNKNOWN, "Rf": UNKNOWN, "Db": UNKNOWN, "Sg": UNKNOWN, "Bh": UNKNOWN,
    "Hs": UNKNOWN, "Mt": UNKNOWN, "Ds": UNKNOWN, "Rg": UNKNOWN, "Cn": UNKNOWN,
    # 7p
    "Nh": UNKNOWN, "Fl": UNKNOWN, "Mc": UNKNOWN, "Lv": UNKNOWN, "Ts": UNKNOWN, "Og": UNKNOWN
}
"""
Wikipedia: "Rayon atomique" - Table 1 (20/02/2025).

Reference :
J. C. Slater, « Atomic Radii in Crystals »,
The Journal of Chemical Physics, vol. 41, no 10, 1964, p. 3199–3204
"""

##### CALCULATED VALUES #####

# Unknown radii are equal to the Bohr radius by default
CALC_UNKNOWN = 53

clementi_et_al_radii_table_pm = {
    # 1s
    "H": 53*NOREF, "He": 31,
    # 2s
    "Li": 167, "Be": 112,
    # 2p
    "B": 87, "C": 67, "N": 56, "O": 48, "F": 42, "Ne": 38,
    # 3s
    "Na": 190, "Mg": 145,
    # 3p
    "Al": 118, "Si": 111, "P": 98, "S": 88, "Cl": 79, "Ar": 71,
    # 4s
    "K": 243, "Ca": 194,
    # 3d
    "Sc": 184, "Ti": 176, "V": 171, "Cr": 166, "Mn": 161,
    "Fe": 156, "Co": 152, "Ni": 149, "Cu": 145, "Zn": 142,
    # 4p
    "Ga": 136, "Ge": 125, "As": 114, "Se": 103, "Br": 94, "Kr": 88,
    # 5s
    "Rb": 265, "Sr": 219,
    # 4d
    "Y": 212, "Zr": 206, "Nb": 198, "Mo": 190, "Tc": 183,
    "Ru": 178, "Rh": 173, "Pd": 169, "Ag": 165, "Cd": 161,
    # 5p
    "In": 156, "Sn": 145, "Sb": 133, "Te": 123, "I": 115, "Xe": 108,
    # 6s
    "Cs": 298, "Ba": 253,
    # 4f
    "La": 226*NOREF, "Ce": 210*NOREF, "Pr": 247, "Nd": 206, "Pm": 205, "Sm": 238, "Eu": 231,
    "Gd": 233, "Tb": 225, "Dy": 228, "Ho": 226, "Er": 226, "Tm": 222, "Yb": 222,
    # 5d
    "Lu": 217, "Hf": 208, "Ta": 200, "W": 193, "Re": 188,
    "Os": 185, "Ir": 180, "Pt": 177, "Au": 174, "Hg": 171,
    # 6p
    "Tl": 156, "Pb": 154, "Bi": 143, "Po": 135, "At": 127, "Rn": 120,
    # 7s
    "Fr": CALC_UNKNOWN, "Ra": CALC_UNKNOWN,
    # 5f
    "Ac": CALC_UNKNOWN, "Th": CALC_UNKNOWN, "Pa": CALC_UNKNOWN, "U": CALC_UNKNOWN, "Np": CALC_UNKNOWN, "Pu": CALC_UNKNOWN, "Am": CALC_UNKNOWN,
    "Cm": CALC_UNKNOWN, "Bk": CALC_UNKNOWN, "Cf": CALC_UNKNOWN, "Es": CALC_UNKNOWN, "Fm": CALC_UNKNOWN, "Md": CALC_UNKNOWN, "No": CALC_UNKNOWN,
    # 6d
    "Lr": CALC_UNKNOWN, "Rf": CALC_UNKNOWN, "Db": CALC_UNKNOWN, "Sg": CALC_UNKNOWN, "Bh": CALC_UNKNOWN,
    "Hs": CALC_UNKNOWN, "Mt": CALC_UNKNOWN, "Ds": CALC_UNKNOWN, "Rg": CALC_UNKNOWN, "Cn": CALC_UNKNOWN,
    # 7p
    "Nh": CALC_UNKNOWN, "Fl": CALC_UNKNOWN, "Mc": CALC_UNKNOWN, "Lv": CALC_UNKNOWN, "Ts": CALC_UNKNOWN, "Og": CALC_UNKNOWN
}
"""
Wikipedia: "Rayons atomiques des éléments (page de données)" - Table 1, column "valeurs calculées" (20/02/2025).

References :
- E. Clementi and D. L. Raimondi,
"Atomic screening constants from SCF functions.",
The Journal of Chemical Physics, vol. 38, no. 11, 1963, p. 2686-2689.

- E. Clementi, D. L. Raimondi and W. P. Reinhardt,
« Atomic Screening Constants from SCF Functions. II. Atoms with 37 to 86 Electrons »,
The Journal of Chemical Physics, vol. 47, no 4, 1967, p. 1300–1307.
"""