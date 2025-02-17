#!/usr/bin/python
"""
Functions to find spacegroup symmetry on pymatgen Structure objects.
"""


import warnings
from typing import List, Union
from functools import partial
from tqdm.contrib.concurrent import process_map

# PYTHON MATERIAL GENOMICS
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmetrizedStructure

# LOCAL IMPORTS
from .common_asserts import check_type, check_num_value
from .redirect import redirect_c_stdout, redirect_c_stderr

########################################
# LOCAL FUNCTIONS

class SymmNotFoundWarning(UserWarning):
    """Warning for symmetry searching returning None."""


########################################

def structure_symmetrizer(
        structure: Structure,
        symprec: float = 0.01,
        angle_tolerance: float = 5.0
    ) -> Union[Structure, SymmetrizedStructure, None]:
    """
    Try to find spacegroup symmetry of a structure using spglib via pymatgen.
    If the first try does not work, it will retry several times with loosened tolerances.

    Parameters:
        structure (Structure):      Pymatgen Structure object to symmetrize.

        valid_tol (float):          Tolerance in relative atomic positions checking
                                    in Angstroms. If the structure contains atoms that
                                    are closer than valid_tol, the function returns None.
                                    If valid_tol = 0.0, distance checking is disabled.
                                    Defaults to 0.0.

        symprec (float):            Initial position tolerance for symmetry detection
                                    in fractional coordinate. Defaults to 0.01, which
                                    works nicely in most cases.
        
        angle_tolerance (float):    Initial angles tolerance for symmetry detection in degrees.
                                    Defaults to 5.0 degrees, which works nicely in most cases.

    Returns:
        A Pymatgen SymmetrizedStructure object if symmetry detection worked properly,
        or the original Structure object if it could not detect any symmetry.
    """
    check_type(structure, "structure", (Structure,))
    check_type(symprec, "symprec", (float,))
    check_num_value(symprec, "symprec", ">=", 0.0)
    check_type(angle_tolerance, "angle_tolerance", (float,))
    check_num_value(angle_tolerance, "angle_tolerance", ">=", 0.0)

    with redirect_c_stdout(None), redirect_c_stderr(None):

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
                warnings.warn(
                        "spglib could not find any symmetry group for this structure, "
                        "maybe it contains some too short interatomic distances.",
                        SymmNotFoundWarning
                )
                return structure

            raise exc

        return sym_struct

########################################

def batch_symmetrizer(
        structures: List[Structure],
        symprec: float = 0.01,
        angle_tolerance: float = 5.0,
        workers: int = 1
    ):
    """
    Use multiprocess to find symmetry spacegroups for a list of pymatgen Structure objects.

    Parameters:
        structures ([Structure]):   A list of all structures that need a symmetry analysis.

        symprec (float):            Initial position tolerance for symmetry detection
                                    in fractional coordinate. Defaults to 0.01, which
                                    works nicely in most cases.
        
        angle_tolerance (float):    Initial angles tolerance for symmetry detection in degrees.
                                    Defaults to 5.0 degrees, which works nicely in most cases.
        
        workers (int):              Number of parallel processes to create.

    Returns:
        A list of either pymatgen SymmetrizedStructure objects
        when symmetry detection worked properly, or the unchanged
        Structure object if symmetry detection failed.
    """
    check_type(structures, "structures", (List,))
    for idx, struct in enumerate(structures):
        check_type(struct, f"structures[{idx}]", (Structure,))
    check_type(symprec, "symprec", (float,))
    check_num_value(symprec, "symprec", ">=", 0.0)
    check_type(angle_tolerance, "angle_tolerance", (float,))
    check_num_value(angle_tolerance, "angle_tolerance", ">=", 0.0)

    nbr_structs = len(structures)
    chunksize   = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)
    set_structure_symmetrizer = partial(
        structure_symmetrizer,
        symprec=symprec,
        angle_tolerance=angle_tolerance
    )

    return list(filter(
        None,
        process_map(
            set_structure_symmetrizer,
            structures,
            max_workers=workers,
            chunksize=chunksize,
            desc="Symmetrize structures"
        )
    ))


########################################
