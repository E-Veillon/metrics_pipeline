#!/usr/bin/python
'''
Global variables for storing experimental fits and values used in the pipeline.
'''


########################################
# NOTE: Unused for now.
E_O2_FIT = -1,36 # eV per O2 in the formation reaction

'''
O2 Energy as fitted by Wang et al.

Reference: 
   L. Wang, T. Maxisch, G. Ceder, Physical Review B 73 (2006) 195107.
'''

########################################

U_VALUES = {
    'F': {
        'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4.0, #'Cu': 4 -> 4.0
        'Fe': 4.0, 'Mn': 3.9, 'Mo': 3.5, 'Nb': 1.5, #'Mo': 4.38 -> 3.5 (according to the reference)
        'Ni': 6.0, 'Re': 2.0, 'Ta': 2.0, 'V': 3.1,  #'Ni': 6 -> 6.0, 'Re': 2 -> 2.0, 'Ta': 2 -> 2.0
        'W': 4.0
    },
    'O': {
        'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4.0, #'Cu': 4 -> 4.0
        'Fe': 4.0, 'Mn': 3.9, 'Mo': 3.5, 'Nb': 1.5, #'Mo': 4.38 -> 3.5 (according to the reference)
        'Ni': 6.0, 'Re': 2.0, 'Ta': 2.0, 'V': 3.1,  #'Ni': 6 -> 6.0, 'Re': 2 -> 2.0, 'Ta': 2 -> 2.0
        'W': 4.0                          
    },
    'S': {
        'Fe': 1.9, 'Mn': 2.5
    }}

'''
Values of the Hubbard U correction used in GGA + U framework, as fitted by Jain et al.

Reference:
    A. Jain, G. Hautier, C.J. Moore, S.P. Ong, 
    C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder, 
    Computational Materials Science, 50, 2295-2310 (2011)
'''

########################################
#NOTE: Unused for now.
DELTA_E_M = {
    'F': {
        'Ag': 0.0, 'Co': 0.0, 'Cr': 0.0, 'Cu': 0.0, 
        'Fe': 0.0, 'Mn': 0.0, 'Mo': 0.0, 'Nb': 0.0, 
        'Ni': 0.0, 'Re': 0.0, 'Ta': 0.0, 'V': 0.0, 
        'W': 0.0
    },
    'O': {
        'Ag': 0.0, 'Co': 0.0, 'Cr': 0.0, 'Cu': 0.0, 
        'Fe': 0.0, 'Mn': 0.0, 'Mo': 0.0, 'Nb': 0.0, 
        'Ni': 0.0, 'Re': 0.0, 'Ta': 0.0, 'V': 0.0, 
        'W': 0.0                          
    },
    'S': {
        'Fe': 0.0, 'Mn': 0.0
}}

'''
Additional correction ΔE_M to consider on GGA + U calculations,
when using the mixed GGA / GGA + U scheme proposed by Jain et al.

Reference:
    A. Jain, G. Hautier, S.P. Ong, C.J. Moore, 
    C.C. Fischer, K.A. Persson, and G. Ceder, 
    Phys. Rev. B, 84, 045115 (2011)
'''

########################################
#NOTE: Unused for now.
EXP_DELTA_H = {
    'F': {
        'Ag': {
            'Ag5F': 0.0, 'Ag2F': 0.0, 'AgF': -1.052, 'Ag2F3': 0.0, 
            'AgF2': 0.0, 'Ag2F5': 0.0, 'Ag3F8': 0.0, 'AgF3': 0.0
        },
        'Co': {
            'Co3F': 0.0, 'Co2F': 0.0, 'CoF': 0.0, 'CoF2': -2.325, 
            'Co2F3': 0.0, 'Co2F5': 0.0, 'CoF3': -2.052, 'CoF4': 0.0, 
            'CoF6': 0.0
        },
        'Cr': {
            'Cr5F': 0.0, 'Cr3F': 0.0, 'Cr2F': 0.0, 'CrF': 0.0, 
            'Cr2F3': 0.0, 'CrF2': -2.701, 'Cr2F5': 0.0, 'CrF3': -3.006, 
            'CrF4': -2.584, 'CrF5': 0.0, 'CrF6': 0.0
        },
        'Cu': {
            'Cu5F': 0.0, 'Cu3F': 0.0, 'Cu2F': 0.0, 'CuF': -1.347, 
            'CuF2': -1.862, 'Cu2F3': 0.0, 'CuF3': 0.0
        },
        'Fe': {
            'Fe5F': 0.0, 'Fe3F': 0.0, 'Fe2F': 0.0, 'Fe3F2': 0.0, 
            'FeF': 0.0, 'FeF2': -2.463, 'Fe2F5': 0.0, 'FeF3': -2.566, 
            'FeF4': 0.0, 'FeF6': 0.0
        },
        'Mn': {
            'Mn8F': 0.0, 'Mn5F': 0.0, 'Mn3F': 0.0, 'Mn2F': 0.0, 
            'MnF': 0.0, 'Mn3F2': 0.0, 'MnF2': -2.954, 'Mn2F5': 0.0, 
            'Mn3F8': 0.0, 'MnF3': -2.775, 'MnF4': -2.276, 'MnF6': 0.0, 
            'MnF7': 0.0
        },
        'Mo': {
            'Mo5F': 0.0, 'Mo3F': 0.0, 'Mo2F': 0.0, 'Mo3F2': 0.0, 
            'MoF': 0.0, 'MoF2': 0.0, 'MoF3': -2.358, 'MoF4': 0.0, 
            'Mo2F9': 0.0, 'MoF5': -2.367, 'MoF6': 0.0
        },
        'Nb': {
            'Nb5F': 0.0, 'Nb3F': 0.0, 'Nb2F': 0.0, 'Nb3F2': 0.0, 
            'NbF': 0.0, 'NbF2': 0.0, 'Nb2F5': 0.0, 'NbF3': 0.0, 
            'NbF4': 0.0, 'NbF5': -3.133
        },
        'Ni': {
            'Ni5F': 0.0, 'Ni3F': 0.0, 'Ni2F': 0.0, 'Ni3F2': 0.0, 
            'NiF': 0.0, 'NiF2': -2.271, 'Ni2F5': 0.0, 'NiF3': 0.0, 
            'NiF4': 0.0, 'NiF6': 0.0
        },
        'Re': {
            'Re5F': 0.0, 'Re3F': 0.0, 'Re2F': 0.0, 'Re3F2': 0.0, 
            'ReF': 0.0, 'ReF2': 0.0, 'ReF3': 0.0, 'ReF4': 0.0, 
            'ReF5': 0.0, 'ReF6': 0.0, 'ReF7': 0.0
        },
        'Ta': {
            'Ta5F': 0.0, 'Ta3F': 0.0, 'Ta2F': 0.0, 'TaF': 0.0, 
            'TaF2': 0.0, 'TaF3': 0.0, 'TaF4': 0.0, 'TaF5': -3.288
        },
        'V': {
            'V5F': 0.0, 'V3F': 0.0, 'V2F': 0.0, 'V3F2': 0.0, 
            'VF': 0.0, 'VF2': -3.428, 'V2F5': 0.0, 'V3F8': 0.0, 
            'VF3': -3.361, 'VF4': -2.909,  'VF5': 0.0
        },
        'W': {
            'W5F': 0.0, 'W3F': 0.0, 'W2F': 0.0, 'W3F2': 0.0, 
            'WF': 0.0, 'WF2': 0.0, 'WF3': 0.0, 'WF4': 0.0, 
            'WF5': 0.0, 'WF6': -2.587
        }
    },
    'O': {
        'Ag': {
            'Ag5O': 0.0, 'Ag3O': 0.0, 'Ag2O': -0.108, 'AgO': 0.0, 
            'Ag4O5': 0.0, 'Ag3O4': 0.0, 'Ag2O3': 0.0, 'AgO2': 0.0, 
            'AgO3': 0.0, 'AgO4': 0.0
        },
        'Co': {
            'Co3O': 0.0, 'Co2O': 0.0, 'CoO': -1.233, 'Co6O7': 0.0, 
            'Co4O5': 0.0, 'Co3O4': -1.356, 'Co2O3': 0.0, 'Co4O7': 0.0, 
            'CoO2': 0.0, 'Co2O5': 0.0, 'Co4O11': 0.0, 'CoO3': 0.0, 
            'CoO4': 0.0
        },
        'Cr': {
            'Cr5O': 0.0, 'Cr3O': 0.0, 'Cr2O': -0.108, 'CrO': 0.0, 
            'Cr4O5': 0.0, 'Cr4O7': 0.0, 'CrO2': -2.009, 'Cr2O5': 0.0, 
            'Cr2O3': -2.364, 'Cr3O8': 0.0, 'Cr8O21': -1.685, 'Cr4O11': 0.0, 
            'CrO3': -1.506, 'CrO4': 0.0
        },
        'Cu': {
            'Cu8O': 0.0, 'Cu5O': 0.0, 'Cu3O': 0.0, 'Cu2O': -0.589, 
            'Cu4O3': 0.0, 'CuO': -0.807, 'Cu4O5': 0.0, 'Cu3O4': 0.0, 
            'Cu2O3': 0.0, 'CuO2': 0.0, 'Cu4O11': 0.0, 'CuO3': 0.0, 
            'CuO4': 0.0
        },
        'Fe': {
            'Fe5O': 0.0, 'Fe3O': 0.0, 'Fe2O': 0.0, 'FeO': -1.415, 
            'Fe4O5': 0.0, 'Fe3O4': -1.660, 'Fe2O3': -1.710, 'Fe5O8': 0.0, 
            'FeO2': 0.0, 'Fe2O5': 0.0, 'Fe4O11': 0.0, 'FeO3': 0.0, 
            'FeO4': 0.0
        },
        'Mn': {
            'Mn8O': 0.0, 'Mn5O': 0.0, 'Mn3O': 0.0, 'Mn2O': 0.0, 
            'MnO': -1.996, 'Mn7O8': 0.0, 'Mn4O5': 0.0, 'M3O4': -2.053, 
            'Mn2O3': -1.983, 'Mn4O7': 0.0, 'Mn5O8': 0.0, 'MnO2': -1.797, 
            'Mn5O12': 0.0, 'Mn2O5': 0.0, 'Mn4O11': 0.0, 'MnO3': 0.0, 
            'MnO4': 0.0
        },
        'Mo': {
            'Mo5O': 0.0, 'Mo3O': 0.0, 'Mo2O': 0.0, 'MoO': 0.0, 
            'Mo4O5': 0.0, 'Mo3O4': 0.0, 'Mo2O3': 0.0, 'MoO2': -2.036, 
            'Mo2O5': 0.0, 'Mo4O11': 0.0, 'Mo9O26': 0.0, 'MoO3': -1.929, 
            'MoO4': 0.0
        },
        'Nb': {
            'Nb3O': 0.0, 'Nb2O': 0.0, 'Nb6O5': 0.0, 'NbO': -2.175, 
            'Nb4O5': 0.0, 'Nb3O4': 0.0, 'Nb2O3': 0.0, 'NB3O5': 0.0, 
            'Nb4O7': 0.0, 'NbO2': -2.746, 'Nb12O29': 0.0, 'Nb2O5': -2.812, 
            'NbO3': 0.0, 'NbO4': 0.0
        },
        'Ni': {
            'Ni5O': 0.0, 'Ni3O': 0.0, 'Ni2O': 0.0, 'Ni3O2': 0.0, 
            'NiO': -1.242, 'Ni7O8': 0.0, 'Ni3O4': 0.0, 'Ni2O3': 0.0, 
            'Ni5O8': 0.0, 'Ni4O7': 0.0, 'NiO2': 0.0, 'Ni4O11': 0.0, 
            'Ni3O8': 0.0, 'NiO3': 0.0, 'NiO4': 0.0
        },
        'Re': {
            'Re5O': 0.0, 'Re3O': 0.0, 'Re2O': 0.0, 'Re3O2': 0.0, 
            'ReO': 0.0, 'Re4O5': 0.0, 'Re2O3': 0.0, 'Re4O7': 0.0, 
            'ReO2': -1.523, 'Re2O5': 0.0, 'ReO3': -1.526, 'Re2O7': -1.455, 
            'ReO4': 0.0
        },
        'Ta': {
            'Ta6O': 0.0, 'Ta5O': 0.0, 'Ta4O': 0.0, 'Ta3O': 0.0, 
            'Ta2O': 0.0, 'Ta3O2': 0.0, 'TaO': 0.0, 'Ta4O5': 0.0, 
            'Ta2O3': 0.0, 'Ta4O7': 0.0, 'TaO2': 0.0, 'Ta2O5': -3.034, 
            'TaO3': 0.0, 'TaO4': 0.0
        },
        'V': {
            'V8O': 0.0, 'V16O3': 0.0, 'V5O': 0.0, 'V3O': 0.0, 
            'V7O3': 0.0, 'V2O': 0.0, 'V3O2': 0.0, 'V5O4': 0.0, 
            'VO': -2.232, 'V7O8': 0.0, 'V4O5': 0.0, 'V3O4': 0.0, 
            'V2O3': -2.522, 'V5O8': 0.0, 'V3O5': 0.0, 'V4O7': 0.0, 
            'V5O9': 0.0, 'V6O11': 0.0, 'V7O13': 0.0, 'V8O15': 0.0, 
            'VO2': -2.475, 'V6O13': 0.0, 'V4O9': 0.0, 'V3O7': 0.0, 
            'V2O5': -2.296, 'V3O8': 0.0, 'V4O11': 0.0, 'VO3': 0.0, 
            'VO4': 0.0
        },
        'W': {
            'W5O': 0.0, 'W3O': 0.0, 'W2O': 0.0, 'W3O2': 0.0, 
            'WO': 0.0, 'W11O12': 0.0, 'W4O5': 0.0, 'W2O3': 0.0, 
            'WO2': -2.038, 'W3O8': 0.0, 'W18O49': 0.0, 'WO3': -2.184, 
            'WO4': 0.0
        }
    },
    'S': {
        'Fe': {
            'Fe8S': 0.0, 'Fe3S': 0.0, 'Fe2S': 0.0, 'Fe5S4': 0.0, 
            'Fe9S8': 0.0, 'FeS': -0.525, 'Fe7S8': 0.0, 'Fe3S4': 0.0, 
            'Fe2S3': 0.0, 'FeS2': -0.588, 'Fe3S8': 0.0, 'FeS3': 0.0, 
            'FeS13': 0.0
        },
        'Mn': {
            'Mn6S': 0.0, 'Mn3S': 0.0, 'Mn2S': 0.0, 'Mn9S8': 0.0, 
            'MnS': -1.110, 'Mn3S4': 0.0, 'Mn2S3': 0.0, 'MnS2': -0.773, 
            'MnS3': 0.0, 'MnS5': 0.0
        }
    }
}

'''
Experimentally measured heats of formations at 298K,
extracted from the Open Quantum Materials Database (OQMD).
These values are used to fit the ΔE_M correction term for GGA + U values 
in the mixed GGA / GGA + U scheme proposed by Jain et al.

NB: A value of 0.0 means the material exists in the database, 
    but no experimental measurement was provided.

References:
  - A. Jain, G. Hautier, S.P. Ong, C.J. Moore, 
    C.C. Fischer, K.A. Persson, and G. Ceder, 
    Phys. Rev. B, 84, 045115 (2011)

  - Saal, J. E., Kirklin, S., Aykol, M., Meredig, B., and Wolverton, C. 
    "Materials Design and Discovery with High-Throughput Density Functional Theory: 
    The Open Quantum Materials Database (OQMD)", JOM 65, 1501-1509 (2013). 
    doi:10.1007/s11837-013-0755-4

  - Kirklin, S., Saal, J.E., Meredig, B., Thompson, A., 
    Doak, J.W., Aykol, M., Rühl, S. and Wolverton, C. 
    "The Open Quantum Materials Database (OQMD): 
    assessing the accuracy of DFT formation energies", 
    npj Computational Materials 1, 15010 (2015). 
    doi:10.1038/npjcompumats.2015.10

OQMD Website: https://www.oqmd.org/
'''

########################################

_EL_PER_XC_VOL_MIN = {
    'LDA_spd': 50, 'PBE_spd': 59, 'AM05_spd': 60, 
    'LDA_sp': 43, 'PBE_sp': 52, 'AM05_sp': 52
}

_EL_PER_XC_VOL_BEST = {
    'LDA_spd': 63, 'PBE_spd': 72, 'AM05_spd': 76, 
    'LDA_sp': 56, 'PBE_sp': 68, 'AM05_sp': 70
}

_EL_PER_XC_VOL_MAX = {
    'LDA_spd': 80, 'PBE_spd': 88, 'AM05_spd': 91, 
    'LDA_sp': 78, 'PBE_sp': 87, 'AM05_sp': 92
}

EL_PER_XC_VOL = {
    'MIN': _EL_PER_XC_VOL_MIN, 
    'BEST': _EL_PER_XC_VOL_BEST, 
    'MAX': _EL_PER_XC_VOL_MAX
}

'''
Values of N* (the number of electrons per exchange-correlation volume) 
used by M.K.Y. Chan and G. Ceder to determine the number 'n' of electrons 
to add to or remove from the simulated structure in the Δ-Sol method.

Reference:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (values in Table I)
'''

########################################
