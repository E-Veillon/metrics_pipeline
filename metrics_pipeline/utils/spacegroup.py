#!/usr/bin/python
"""
Functions to find spacegroup symmetry on pymatgen Structure objects.
"""


import warnings
from typing import List, Union
from functools import partial
from tqdm import tqdm
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

SPG_NUM_TO_PG = {
    # Triclinic
    1: "1", 2: "-1",
    # Monoclinic
    **dict.fromkeys(list(range(3, 6)), "2"),
    **dict.fromkeys(list(range(6, 10)), "m"),
    **dict.fromkeys(list(range(10, 16)), "2/m"),
    # Orthorhombic
    **dict.fromkeys(list(range(16, 25)), "222"),
    **dict.fromkeys(list(range(25, 47)), "mm2"),
    **dict.fromkeys(list(range(47, 75)), "mmm"),
    # Tetragonal
    **dict.fromkeys(list(range(75, 81)), "4"),
    **dict.fromkeys(list(range(81, 83)), "-4"),
    **dict.fromkeys(list(range(83, 89)), "4/m"),
    **dict.fromkeys(list(range(89, 99)), "422"),
    **dict.fromkeys(list(range(99, 111)), "4mm"),
    **dict.fromkeys(list(range(111, 123)), "-42m"),
    **dict.fromkeys(list(range(123, 143)), "4/mmm"),
    # Trigonal
    **dict.fromkeys(list(range(143, 147)), "3"),
    **dict.fromkeys(list(range(147, 149)), "-3"),
    **dict.fromkeys(list(range(149, 156)), "32"),
    **dict.fromkeys(list(range(156, 162)), "3m"),
    **dict.fromkeys(list(range(162, 168)), "-3m"),
    # Hexagonal
    **dict.fromkeys(list(range(168, 174)), "6"),
    **dict.fromkeys(list(range(174, 175)), "-6"),
    **dict.fromkeys(list(range(175, 177)), "6/m"),
    **dict.fromkeys(list(range(177, 183)), "622"),
    **dict.fromkeys(list(range(183, 187)), "6mm"),
    **dict.fromkeys(list(range(187, 191)), "-6m2"),
    **dict.fromkeys(list(range(191, 195)), "6/mmm"),
    # Cubic
    **dict.fromkeys(list(range(195, 200)), "23"),
    **dict.fromkeys(list(range(200, 207)), "m-3"),
    **dict.fromkeys(list(range(207, 215)), "432"),
    **dict.fromkeys(list(range(215, 221)), "-43m"),
    **dict.fromkeys(list(range(221, 231)), "m-3m"),
}
PG_TO_SYSTEM = {
    **dict.fromkeys(("1", "-1"), "triclinic"),
    **dict.fromkeys(("2", "m", "2/m"), "monoclinic"),
    **dict.fromkeys(("222", "mm2", "mmm"), "orthorhombic"),
    **dict.fromkeys(("4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm"), "tetragonal"),
    **dict.fromkeys(("3", "-3", "32", "3m", "-3m"), "trigonal"),
    **dict.fromkeys(("6", "-6", "6/m", "622", "6mm", "-62m", "6/mmm"), "hexagonal"),
    **dict.fromkeys(("23", "m-3", "432", "-43m", "m-3m"), "cubic"),
}

########################################

def structure_symmetrizer(
        structure: Structure,
        symprec: float = 0.01,
        angle_tolerance: float = 5.0
    ) -> Union[Structure, SymmetrizedStructure]:
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
        workers: int|None = None,
        sequential: bool = False
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

        sequential (bool):          Whether to use sequential for-loop instead of multiprocessing
                                    scheme. If set to True, the 'workers' arg is ignored.
                                    Defaults to False.

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
    if workers is not None:
        check_type(workers, "workers", (int,))
        check_num_value(workers, "workers", ">", 0)

    nbr_structs = len(structures)
    chunksize   = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)
    set_structure_symmetrizer = partial(
        structure_symmetrizer,
        symprec=symprec,
        angle_tolerance=angle_tolerance
    )

    if sequential:
        sym_structs = []
        for struct in tqdm(structures, desc="Symmetrize structures"):
            sym_structs.append(set_structure_symmetrizer(struct))

    else:
        sym_structs = list(filter(
            None,
            process_map(
                set_structure_symmetrizer,
                structures,
                max_workers=workers,
                chunksize=chunksize,
                desc="Symmetrize structures"
            )
        ))
    return sym_structs


########################################
