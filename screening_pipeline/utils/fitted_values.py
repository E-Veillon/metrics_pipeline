'''
Global variables for storing experimental fits and values used in the pipeline.
'''


########################################

E_O2_FIT = 1 #TODO: Modify it when the calculation is done

'''
O2 Energy as fitted with Wang et al. method.

Reference: 
   L. Wang, T. Maxisch, G. Ceder, 
   Physical Review B 73 (2006) 195107.
   (reference 14.57 in screening_pipeline/Bibliography)
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
Values of the Hubbard U correction in GGA + U framework.

Reference:
    A. Jain, G. Hautier, C.J. Moore, S.P. Ong, 
    C.C. Fischer, T. Mueller, K.A. Persson, and G. Ceder, 
    Computational Materials Science, 50, 2295-2310 (2011)
    (reference 14 in screening_pipeline/Bibliography)
'''

########################################

DELTA_E_M = { #TODO: Modify it when values are fitted
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
    (reference 31 in screening_pipeline/Bibliography)
'''

########################################

EXP_DELTA_H = { #TODO: Modify it when values are found
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
            'TaF2': 0.0, 'TaF3': 0.0, 'TaF4': 0.0, 'TaF4': 0.0, 
            'TaF5': -3.288
        }, 
        'V': {
            'V5F': 0.0, 'V3F': 0.0, 'V2F': 0.0, 'V3F2': 0.0, 
            'VF': 0.0, 'VF2': -3.428, 'V2F5': 0.0, 'V3F8': 0.0, 
            'VF3': -3.361, 'VF4': -2.909, 'VF4': 0.0, 'VF5': 0.0
        }, 
        'W': {
            'W5F': 0.0, 'W3F': 0.0, 'W2F': 0.0, 'W3F2': 0.0, 
            'WF': 0.0, 'WF2': 0.0, 'WF3': 0.0, 'WF4': 0.0, 
            'WF4': 0.0, 'WF5': 0.0, 'WF6': -2.587
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
        'Cr': {}, 
        'Cu': {}, 
        'Fe': {'FeO': -1.415, 'Fe3O4': -1.660, 'Fe2O3': -1.710}, 
        'Mn': {}, 
        'Mo': {}, 
        'Nb': {}, 
        'Ni': {}, 
        'Re': {}, 
        'Ta': {}, 
        'V': {}, 
        'W': {}                          
    }, 
    'S': {
        'Fe': {}, 
        'Mn': {}
    }}

'''
Experimentally measured heats of formations at 298K,
extracted from the Open Quantum Materials Database (OQMD).
These values are used to fit the ΔE_M correction term for GGA + U values 
in the mixed GGA / GGA + U scheme proposed by Jain et al.

NB: A value of 0.0 means the material exists in the database, 
    but no experimental measurement was done yet.

References:
  - A. Jain, G. Hautier, S.P. Ong, C.J. Moore, 
    C.C. Fischer, K.A. Persson, and G. Ceder, 
    Phys. Rev. B, 84, 045115 (2011)
    (reference 31 in screening_pipeline/Bibliography)

  - Saal, J. E., Kirklin, S., Aykol, M., Meredig, B., and Wolverton, C. 
    "Materials Design and Discovery with High-Throughput Density Functional Theory: 
    The Open Quantum Materials Database (OQMD)", JOM 65, 1501-1509 (2013). 
    doi:10.1007/s11837-013-0755-4
    (reference OQMD1 in screening_pipeline/Bibliography)

  - Kirklin, S., Saal, J.E., Meredig, B., Thompson, A., 
    Doak, J.W., Aykol, M., Rühl, S. and Wolverton, C. 
    "The Open Quantum Materials Database (OQMD): 
    assessing the accuracy of DFT formation energies", 
    npj Computational Materials 1, 15010 (2015). 
    doi:10.1038/npjcompumats.2015.10
    (reference OQMD2 in screening_pipeline/Bibliography)

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
used by M.K.Y. Chan and G. Ceder in their presentation of Δ-Sol method.

Reference:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (reference 32 in screening_pipeline/Bibliography, values in Table I)
'''

########################################