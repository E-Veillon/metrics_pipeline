'''Functions used in the script perturbator.py.'''


########################################
# SYSTEM I/O MODULES

from typing import Tuple, List, Union, Dict, Sequence, Any
import re

########################################
# OPTIMIZATION MODULES

import itertools
from functools import partial
from random import uniform, choice
from tqdm.contrib.concurrent import process_map

########################################
# PYTHON MATERIALS GENOMICS PACKAGE

from pymatgen.core import Structure, PeriodicSite, Lattice
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SymmetrizedStructure

########################################
# CIF I/O FUNCTIONS

def extract_cif_from_file(filename: str) -> List[str]:
    """
    Separates concatenated structures from a cif file.

    Parameters:
        filename (str):         Path to a cif file.

        keep_rare_gases (bool): Whether the structures containing rare gases should be kept or not. 
                                Defaults to false.
    
    Returns:
        List[str]: List of CIF strings, each one representing a single structure.
    """

    # load file
    with open(filename, "r") as input_file:
        data = input_file.read()

    # match structures data with regular expression
    matcher = re.compile(r"^data.*?$(?=\ndata|\Z)", re.MULTILINE | re.DOTALL)
    cif_str_structs = matcher.findall(data)

    return cif_str_structs


def cif_str_to_struct(cif_str: str) -> Structure:
    """
    Parses data from a cif formatted string and converts it to a pymatgen Structure object.

    Parameters:
        cif_str (str): The string of a structure encoded in the cif format.

    Returns:
        A Pymatgen Structure object.
    """

    parsed_str  = CifParser.from_str(cif_string=cif_str)
    struct_list = parsed_str.parse_structures(primitive=False)
    structure   = struct_list[0]
    return structure


def struct_to_cif_str(
        structure: Union[Structure, SymmetrizedStructure], 
        significant_figures: int = 8
    ) -> str:
    '''
    Converts a pymatgen Structure object into a CIF formatted string.

    Parameters:
        structure (Structure): Structure object to convert.
    
    Returns:
        str: CIF formatted string.
    '''

    if not isinstance(structure, Structure):
        raise TypeError('Cannot write CIF data for a non-structure object.')

    cif_writer = CifWriter(
        struct=structure, 
        significant_figures=significant_figures
    )

    if not isinstance(structure, SymmetrizedStructure):
        cif_str = str(cif_writer)
        return cif_str
    
    # Extract the data dict from CifWriter
    cif_block_list = list(cif_writer.cif_file.data.values())
    cif_block      = cif_block_list[0]
    data_dict      = cif_block.data

    # Extract symmetry infos from the SymmetrizedStructure object
    struct_spg     = structure.spacegroup
    xyz_ops        = [op.as_xyz_str() for op in struct_spg]
    equiv_sites    = structure.equivalent_sites

    # Get the ordering for species listing

    # Inside a list of equivalent sites, we only keep the one with minimal fractional coordinates, 
    # and take the number of equivalent sites as the multiplicity
    frac_coord_sorting_key = lambda s: tuple(abs(x) for x in s.frac_coords)
    unique_sites: List[Tuple[PeriodicSite, int]] = [(sorted(sites, key=frac_coord_sorting_key)[0],len(sites)) for sites in equiv_sites]
    
    # Between non-equivalent sites, we sort them firstly by ascending electronegativity, 
    # then by descending multiplicity, then by ascending frac coordinates
    electroneg_sorting_key = lambda t: (t[0].species.average_electroneg,-t[1],t[0].a,t[0].b,t[0].c)
    sorted_unique_sites: List[Tuple[PeriodicSite, int]] = sorted(unique_sites, key=electroneg_sorting_key)

    atom_site_type_symbol   = []
    atom_site_symmetry_mult = []
    atom_site_fract_x       = []
    atom_site_fract_y       = []
    atom_site_fract_z       = []
    atom_site_label         = []
    atom_site_occupancy     = []

    format_str = f"{{:.{significant_figures}f}}"
    count      = 0

    for site, mult in sorted_unique_sites:
        for specie, occupancy in site.species.items():
            atom_site_type_symbol.append(str(specie))
            atom_site_symmetry_mult.append(f"{mult}")
            atom_site_fract_x.append(format_str.format(site.a))
            atom_site_fract_y.append(format_str.format(site.b))
            atom_site_fract_z.append(format_str.format(site.c))
            site_label = site.label if site.label != site.species_string else f"{specie.symbol}{count}"
            atom_site_label.append(site_label)
            atom_site_occupancy.append(str(occupancy))
            count += 1

    data_dict['_symmetry_space_group_name_H-M']   = struct_spg.int_symbol
    data_dict['_symmetry_Int_Tables_number']      = struct_spg.int_number
    data_dict['_symmetry_equiv_pos_site_id']      = [f'{idx}' for idx in range(1, len(xyz_ops) + 1)]
    data_dict['_symmetry_equiv_pos_as_xyz']       = xyz_ops
    data_dict["_atom_site_type_symbol"]           = atom_site_type_symbol
    data_dict["_atom_site_label"]                 = atom_site_label
    data_dict["_atom_site_symmetry_multiplicity"] = atom_site_symmetry_mult
    data_dict["_atom_site_fract_x"]               = atom_site_fract_x
    data_dict["_atom_site_fract_y"]               = atom_site_fract_y
    data_dict["_atom_site_fract_z"]               = atom_site_fract_z
    data_dict["_atom_site_occupancy"]             = atom_site_occupancy

    cif_str = str(cif_block)

    return cif_str

#def read_cif(
#    filename: str,
#    workers: int = 1,
#    keep_rare_gases: bool = False, 
#    keep_rare_earths: bool = False
#) -> Tuple[List[Structure], int, int]:
#    """
#    Reads a cif file containing concatenated structures data and decode them using multiprocess.
#
#    Parameters:
#        filename (str):             Name of the input CIF file.
#        
#        workers (int):              Number of processes to use in parallel.
#        
#        keep_rare_gases (bool):     Whether structures containing rare gases should be kept.
#                                    Defaults to false.
#        
#        keep_rare_earths (bool):    Whether structures containing f-block elements should be kept.
#                                    Defaults to False.
#    
#    Returns:
#        List[Structure]: Decoded Structure objects in a list.
#        Int: Number of structures containing rare gases discarded.
#        Int: Number of structures containing rare earth elements discarded.
#    """
#
#    struct_strings = extract_cif_from_file(filename)
#    nbr_rare_gas_structs, nbr_rare_earth_structs = 0, 0
#
#    if not keep_rare_gases:
#        struct_strings, nbr_rare_gas_structs = discard_rare_gas_structures(struct_strings)
#        print(f'{nbr_rare_gas_structs} rare gas structures ignored')
#    
#    if not keep_rare_earths:
#        struct_strings, nbr_rare_earth_structs = discard_rare_earth_structures(struct_strings)
#        print(f'{nbr_rare_earth_structs} rare earths structures ignored')
#
#    nbr_struct = len(struct_strings)
#    assert nbr_struct > 0, \
#    "No structure data found in provided file (maybe they all have been discarded ?)"
#
#    chunksize = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)
#
#    structs_list = list(process_map(
#        cif_str_to_struct,
#        struct_strings,
#        max_workers=workers,
#        chunksize=chunksize,
#        desc="load and read data",
#    ))
#
#    return structs_list, nbr_rare_gas_structs, nbr_rare_earth_structs


def simply_read_cif(
    filename: str,
    workers: int = 1
) -> Tuple[List[Structure], int, int]:
    """
    Reads a cif file containing concatenated structures data and decode them using multiprocess.
    Simplified version of read_cif() that cannot do anything else than Structure conversion.

    Parameters:
        filename (str):             Name of the input CIF file.
        
        workers (int):              Number of processes to use in parallel.

    Returns:
        List[Structure]: Decoded Structure objects in a list.
    """

    struct_strings = extract_cif_from_file(filename)

    nbr_struct = len(struct_strings)

    assert nbr_struct > 0, "No structure data found in provided file."

    chunksize = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    structs_list = list(
        process_map(
            cif_str_to_struct,
            struct_strings,
            max_workers=workers,
            chunksize=chunksize,
            desc="load and read data",
        )
    )
    return structs_list


def write_cif(
    filename: str,
    structures: List[Structure],
    workers: int = 1,
) -> None:
    """
    Encode multiple structures in CIF format and write them in a file using multiprocess.

    Parameters:
        filename (str):               Name of the input file.

        structures (List[Structure]): The structures to encode.

        workers (int):                Number of parallel processes to use.
    """

    nbr_struct = len(structures)
    chunksize  = (min(nbr_struct // 100, 10) if nbr_struct >= 200 else 1)

    encoded_cif = list(
        process_map(
            struct_to_cif_str, 
            structures,
            max_workers=workers, 
            chunksize=chunksize, 
            desc="Writing data in cif format"
        )
    )

    with open(filename, "wt") as out_file:
        if all([cif_str.endswith("\n") for cif_str in encoded_cif]):
            out_file.write("".join(encoded_cif))
        else:
            out_file.write("\n".join(encoded_cif))


########################################
# OTHER FUNCTIONS


def flatten(sequence: Sequence[Any], level_of_flattening: int = 1) -> List:
    '''
    Unpacks a nested sequence without modifying elements order.

    Parameters:
        sequence (Sequence):        The sequence to unpack.

        level_of_flattening (Int):  The number of nested levels to unpack.
                                    Defaults to 1.
    
    Returns:
        A list flattened the specified number of times.
    '''

    for _ in range(1, level_of_flattening + 1):
        sequence = list(itertools.chain.from_iterable(sequence))
    return sequence


########################################


def get_random_sign() -> int:
    """Returns randomly +1 or -1."""
    return choice((-1, 1))


########################################


def get_perturbs_dict(
        parameters: List[str],
        min_perturbs: List[float],
        max_perturbs: List[float]
    ) -> Dict[str, Union[Dict[str, float], None]]:
    """
    Parse parameters arguments and values into a dict compiling all perturbations to apply.
    
    Parameters:
        parameters ([str]): Which structure parameters have to be perturbed.
                            Supported parameters:
                            - 'sites' for site positions,
                            - 'lengths' for all lattice vectors lengths (a, b, c),
                            - 'angles' for all lattice angles (alpha, beta, gamma),
                            - 'a', 'b', 'c', 'alpha', 'beta', 'gamma' for individual lattice parameters.
        
        min_perturbs ([float]): Minimal possible perturbation for the parameter at same index in 'parameters'.

        max_perturbs ([float]): Maximal possible perturbation for the parameter at same index in 'parameters'.
    
    Returns:
        A dict with a key for each individual parameter containing corresponding min and max possible perturbation.
    """

    perturbs_dict = dict.fromkeys(("sites", "a", "b", "c", "alpha", "beta", "gamma"))

    for param, mini, maxi in zip(parameters, min_perturbs, max_perturbs):
        match param:
            case "lengths":
                perturbs_dict.update({
                    "a": {"min": mini, "max": maxi},
                    "b": {"min": mini, "max": maxi},
                    "c": {"min": mini, "max": maxi}
                })
            case "angles":
                perturbs_dict.update({
                    "alpha": {"min": mini, "max": maxi},
                    "beta": {"min": mini, "max": maxi},
                    "gamma": {"min": mini, "max": maxi}
                })
            case "sites"|"a"|"b"|"c"|"alpha"|"beta"|"gamma":
                perturbs_dict.update({param: {"min": mini, "max": maxi}})
            case str(): raise ValueError(f"Unsupported parameter '{param}'.")
            case _: raise TypeError(f"Wrong type in for parameter '{param}'. Expected 'str', got '{type(param)}'.")

    return perturbs_dict


########################################


def _perturb_site_positions(structure: Structure, min: float, max: float) -> Structure:
    """
    Perturb randomly all site positions of a structure.
    Perturbations are ranging between values given in 'min' and 'max' arguments, in angstroms.
    """
    perturbed_struct = structure.copy()
    perturbed_struct.perturb(min_distance=min,distance=max)
    return perturbed_struct


########################################


def is_valid_lattice(lattice: Lattice) -> bool:
    """Check whether given lattice has physically valid parameters values."""
    a0 = 0.53 # a0 = Bohr radius, approximately equal to 0.53 angstroms.
    min_vol = (2 * a0)**3 # Minimal cell volume of (2 * a0)^3 angstroms^3.
    min_len = 2 * a0 # a, b, c must be superior than 2 * a0 angstroms.
    min_angle = 5.0 # alpha, beta, gamma must be superior than 5.0 degrees.

    return (
        lattice.volume >= min_vol
        and all([length >= min_len for length in lattice.lengths]) 
        and all([angle >= min_angle for angle in lattice.angles])
    )


########################################


def _perturb_lattice_parameters(
        structure: Structure, perturbs_dict: Dict, retries: int = 4) -> Structure:
    """
    Perturb randomly lattice parameters of a structure.
    Perturbations are ranging between values given in 'perturbs_dict'.
    Perturbations are expressed in angstroms for lengths and degrees for angles.
    """
    old_lattice_params = structure.lattice.params_dict
    new_lattice_params = dict.fromkeys(("a", "b", "c", "alpha", "beta", "gamma"))

    for param in old_lattice_params:
        if perturbs_dict.get(param) is None:
            new_lattice_params[param] = old_lattice_params[param]
        else:
            perturb = get_random_sign() * uniform(
                perturbs_dict[param].get("min"),
                perturbs_dict[param].get("max")
            )
            new_lattice_params[param] = old_lattice_params[param] + perturb

    new_lattice = Lattice.from_parameters(
        a=new_lattice_params["a"],
        b=new_lattice_params["b"],
        c=new_lattice_params["c"],
        alpha=new_lattice_params["alpha"],
        beta=new_lattice_params["beta"],
        gamma=new_lattice_params["gamma"]
    )
    # As the process is random, it might sometimes return an unphysical lattice.
    # If it is the case, try to get another random perturbation which leads to
    # physically valid values.
    if not is_valid_lattice(new_lattice):
        if retries > 0:
            print(f"Unphysical lattice generated, retrying... (retries left: {retries})")
            retries -= 1
            return _perturb_lattice_parameters(structure, perturbs_dict, retries=retries)
        else:
            raise ValueError(
                "The following structure's lattice could not be perturbed "
                "in a physically viable way, you should check if the original "
                "structure is physical and/or if given lattice perturbation amplitudes "
                "makes it possible to have a physical lattice "
                "(i.e. volume > 1.2 Angstroms^3; a, b and c > 1.06 Angstroms; "
                "angles > 1.0 degree):\n"
                f"{structure}\n"
                "---Given lattice perturbation parameters---\n"
                f"on 'a' length: {perturbs_dict['a']}\n"
                f"on 'b' length: {perturbs_dict['b']}\n"
                f"on 'c' length: {perturbs_dict['c']}\n"
                f"on 'alpha' angle: {perturbs_dict['alpha']}\n"
                f"on 'beta' angle: {perturbs_dict['beta']}\n"
                f"on 'gamma' angle: {perturbs_dict['gamma']}\n"
            )

    perturbed_struct = structure.copy()
    perturbed_struct.lattice = new_lattice

    return perturbed_struct


########################################


def generate_perturbed_structs(
        structure: Structure,
        perturbs_dict: Dict[str, Union[Dict[str, float], None]], 
        sample_size: int = 1,
        modified_lattice: bool = False,
        lattice_retries: int = 4
    ) -> List[Structure]:
    """
    Generate perturbed structures from a reference.

    Parameters:
        structure (Structure):      A reference structure to perturb.

        perturbs_dict (Dict):       A dict containing min and max boundaries of
                                    perturbation for each structure parameter
                                    ("sites", "a", "b", "c", "alpha", "beta", "gamma").
        
        sample_size (int):          Number of perturbed structures should be
                                    generated for each input structure. Defaults to 1.
        
        modified_lattice (bool):    Whether lattice parameters should be perturbed.
                                    Defaults to False.

        lattice_retries (int):      Number of times the lattice perturbations
                                    can be retried when they result in an unphysical
                                    lattice before raising an error. Defaults to 4.
        
        Returns:
            The list of generated structures.
    """
    assert isinstance(structure, Structure), (
        f"'structure' argument value must be a Structure object, got '{type(structure)}' instead."
    )
    assert isinstance(perturbs_dict, Dict), (
        f"'perturbs_dict' argument value must be a dict, got '{type(perturbs_dict)}' instead."
    )
    assert isinstance(sample_size, int), (
        f"'sample_size' argument value must be a int, got '{type(sample_size)}' instead."
    )
    assert sample_size >= 1, "'sample_size' argument value must be strictly positive."

    perturbed_structs = []

    for _ in range(sample_size):

        perturbed_struct = None

        # Gestion des positions des sites
        if perturbs_dict.get("sites") is not None:
            perturbed_struct = _perturb_site_positions(
                structure=structure,
                min=perturbs_dict["sites"].get("min"),
                max=perturbs_dict["sites"].get("max")
            )

        # Gestion des paramètres de maille
        if modified_lattice and perturbed_struct is not None:
            perturbed_struct = _perturb_lattice_parameters(
                structure=perturbed_struct,
                perturbs_dict=perturbs_dict,
                retries=lattice_retries
            )
        elif modified_lattice and perturbed_struct is None:
            perturbed_struct = _perturb_lattice_parameters(
                structure=structure,
                perturbs_dict=perturbs_dict,
                retries=lattice_retries
            )

        # Addition du résultat à la liste
        if perturbed_struct is None:
            perturbed_structs.append(structure)

        else:
            perturbed_structs.append(perturbed_struct)

    return perturbed_structs


########################################


def batch_generate_perturbed_structs(
        structures: Sequence[Structure],
        perturbs_dict: Dict[str, Union[Dict[str, float], None]],
        sample_size: int = 1,
        modified_lattice: bool = False,
        lattice_retries: int = 4,
        workers: int = 1
    ) -> List[Structure]:
    """
    Parallelized version of generate_perturbed_structs() designed for a batch of structures.
    
    Parameters:
        structures ([Structure]):   A sequence of structures to perturb.

        perturbs_dict (Dict):       A dict containing min and max boundaries of
                                    perturbation for each structure parameter
                                    ("sites", "a", "b", "c", "alpha", "beta", "gamma").
        
        sample_size (int):          Number of perturbed structures that should be
                                    generated for each input structure. Defaults to 1.
        
        modified_lattice (bool):    Whether lattice parameters should be perturbed.
                                    Defaults to False.

        lattice_retries (int):      Number of times the lattice perturbations
                                    can be retried when they result in an unphysical
                                    lattice before raising an error. Defaults to 4.
        
        workers (int):              Number of parallel processes to spawn.
                                    Defaults to 1.

        Returns:
            The list of all generated structures.
    """
    assert all([isinstance(struct, Structure) for struct in structures]), (
        "Some of the given objects in 'structure' argument are not Structure objects.\n"
        "See below all encountered object types:\n"
        f"{set([type(struct) for struct in structures])}"
    )
    assert isinstance(perturbs_dict, Dict), "'perturbs_dict' argument value must be a dict."
    assert workers >= 1, "'workers' argument value must be strictly positive."

    nbr_structs = len(structures)
    chunksize   = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)
    perturb_setup = partial(generate_perturbed_structs,
        perturbs_dict=perturbs_dict, sample_size=sample_size, 
        modified_lattice=modified_lattice, lattice_retries=lattice_retries
    )

    perturbed_structs = flatten(
        list(
            process_map(
                perturb_setup,
                structures,
                max_workers=workers,
                chunksize=chunksize,
                desc=f"generating {sample_size} perturbed struct(s) per structure"
            )
        )
    )

    return perturbed_structs


########################################