'''
Functions to find spacegroup symmetry on pymatgen Structure objects.
'''


########################################
# TYPE HINTING

from typing import Tuple, List

########################################
# OPTIMIZATION MODULES

from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmetrizedStructure

########################################


def get_default_symmetry(structure: Structure) -> Structure:
    #TODO: modify CifBlock data dict to set it to P1 and spg number 1
    return structure

def _retry_get_symmetrized_structure(
        structure: Structure,
        symprec: float = 0.01,
        angle_tolerance: float = 5.0,
    ) -> Structure:
    '''
    Retry function for symmetry detection with several sets of loosened tolerances,
    called in case it did not work with initial given tolerances.

    Parameters:
        structure (Structure):      Pymatgen Structure object to symmetrize.

        symprec (float):            Initial position tolerance for symmetry detection in fractional coordinate.
                                    Defaults to 0.01, which works nicely in most cases.
        
        angle_tolerance (float):    Initial angles tolerance for symmetry detection in degrees.
                                    Defaults to 5.0 degrees, which works nicely in most cases.

    Returns:
        A Pymatgen SymmetrizedStructure object if symmetry detection worked properly,
        or a Structure object with default P1 spacegroup if it could not detect any.
    '''
    for precision_factor in [2, 3, 5, 10]:

        symmetrizer = SpacegroupAnalyzer(
            structure=structure, 
            symprec=precision_factor*symprec, 
            angle_tolerance=precision_factor*angle_tolerance
        )

        try:
            sym_struct = symmetrizer.get_symmetrized_structure()

        except TypeError:

            if symmetrizer.get_symmetry_dataset() is not None:
                return structure
            
        else:
            return sym_struct
    
    return get_default_symmetry(structure)

def structure_symmetrizer(
        structure: Structure, 
        symprec: float = 0.01, 
        angle_tolerance: float = 5.0
    ) -> Structure | SymmetrizedStructure:
    '''
    Try to find spacegroup symmetry of a structure using spglib via pymatgen.
    If the first try does not work, it will retry several times with loosened tolerances.

    Parameters:
        structure (Structure):      Pymatgen Structure object to symmetrize.

        symprec (float):            Initial position tolerance for symmetry detection in fractional coordinate.
                                    Defaults to 0.01, which works nicely in most cases.
        
        angle_tolerance (float):    Initial angles tolerance for symmetry detection in degrees.
                                    Defaults to 5.0 degrees, which works nicely in most cases.

    Returns:
        A Pymatgen SymmetrizedStructure object if symmetry detection worked properly,
        or a Structure object with default P1 spacegroup if it could not detect any.
    '''
    struct_analyzer = SpacegroupAnalyzer(
        structure=structure, 
        symprec=symprec, 
        angle_tolerance=angle_tolerance
    )
    try:
        sym_struct = struct_analyzer.get_symmetrized_structure()

    except TypeError as exc:
        spglib_result = struct_analyzer.get_symmetry_dataset()

        if spglib_result is None:

            return _retry_get_symmetrized_structure(
                structure=structure, 
                symprec=symprec, 
                angle_tolerance=angle_tolerance
            )

        raise exc

    return sym_struct

def batch_symmetrizer(
        structures: List[Structure], 
        symprec: float = 0.01, 
        angle_tolerance: float = 5.0, 
        workers: int = 1
    ):
    '''
    Use multiprocess to find symmetry spacegroups for a list of pymatgen Structure objects.

    Parameters:
        structures (List[Structure]):   A list of all structures that need a symmetry analysis.

        symprec (float):                Initial position tolerance for symmetry detection in fractional coordinate.
                                        Defaults to 0.01, which works nicely in most cases.
        
        angle_tolerance (float):        Initial angles tolerance for symmetry detection in degrees.
                                        Defaults to 5.0 degrees, which works nicely in most cases.
        
        workers (int):                  Number of parallel processes to create.

    Returns:
        A list containing either pymatgen SymmetrizedStructure objects when symmetry detection worked properly,
        or Structure objects with default P1 spacegroup if detection could not detect any symmetry.
    '''
    def feed_args(structures, symprec, angle_tolerance) -> List[Tuple]:
        return [(struct, symprec, angle_tolerance) for struct in structures]

    nbr_structs = len(structures)
    chunksize   = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)
    
    return list(
        process_map(
            structure_symmetrizer, 
            feed_args(structures, symprec, angle_tolerance), 
            max_workers=workers, 
            chunksize=chunksize
    ))