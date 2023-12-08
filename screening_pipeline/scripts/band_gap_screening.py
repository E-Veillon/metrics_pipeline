#!/usr/bin/python
'''
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (reference 32 in screening_pipeline/Bibliography)
'''

########################################
# SYSTEM I/O MODULES

from datetime import datetime
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter
from pathlib import Path

########################################
# PYTHON MATERIALS GENOMICS PACKAGE

#from pymatgen.io.vasp.outputs import Oszicar, Chgcar

########################################
# LOCAL MODULES

from screening_pipeline.utils.cif_io import read_cif
from screening_pipeline.utils.periodic_table import get_all_valence_electrons, get_delta_sol_el_ratio
from screening_pipeline.utils.vasp_io import vasp_input_files_settings, vasp_launcher

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:
    
    assert Path(args.input_dir).is_dir(), \
    f'{args.input_dir}: No directory found.'

    assert Path(args.output).is_dir(), \
    f'{args.output}: No directory found.'

    assert args.workers >= 1, \
    'The number of workers cannot be negative or zero.'

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK
    
    prog_name = 'band_gap_screening'
    prog_desc = '''
        A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

        Reference for Δ-Sol method:
            M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
            (reference 32 in screening_pipeline/Bibliography)
        '''
    prog_missing_steps = '''
        Missing steps to complete this script:
            - Add necessary in-line command arguments
            - Calculate N0, number of valence electrons in unit cell - OK
            - Choose right N*best (PBE[sp] or PBE[spd]) - OK
            - n = N0/N*best - OK
            - Calculate E(N0), E(N0 + n), E(N0 - n) with VASP (static ionically)
            - Calculate E(gap) = [E(N0 + n) + E(N0 - n) - 2*E(N0)]/n
            - For photocatalysts article: keep materials with 1.3 < E(gap) < 3.6 eV by default
            - Set an option to enable incertainty calculation on the band gap using N*min and N*max
            - This option should let user choose if they want to keep materials according to the uncertainty case:
                * If E(gap) is inside the goal but uncertainty gets out ?
                * If E(gap) is outside the goal but uncertainty gets in ?
                Possible options: keep, keep_aside, discard, auto (keep_aside if at least half the interval is in, else discard)
        '''
    helper_format = RawDescriptionHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_desc, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        'input_dir',
        type=str,
        help='Base directory containing structure directories.', 
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default='./',
        help='''Path to the output directory where VASP files will be written.
                A subdirectory will be created in output directory
                for each structure processed.''', 
        metavar='outdir'
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to create',
        metavar='int',
    )
    args: Namespace = parser.parse_args()

    assert_args(args)
    
    input_dir  = Path(args.input_dir)
    outdir     = Path(args.output)
    workers    = args.workers


    # MAIN BLOCK

    # Créer une fonction qui va chercher les fichiers OSZICAR dans input_path, 
    # puis qui les convertit en objets OSZICAR dont on peut tirer l'énergie E(N0).
    # Attention, il faut garder l'ordre des structures à l'aide d'un dict {num_struct: E(N0)}.
    from typing import Union, Dict, Tuple
    from functools import partial
    from tqdm.contrib.concurrent import process_map
    from pymatgen.core.structure import Structure
    from pymatgen.io.vasp import Poscar, Chgcar, Oszicar
    PathLike = Union[str, Path]

    def extract_vasp_data(
            struct_dir: PathLike = '.', 
            ignore_file: str = 'rejected.txt'
    ) -> Tuple[str, Dict[Structure, Chgcar, float]]:
        '''
        Extracts VASP data from a previous run for one structure. 
        Keeps only relevant data for Δ-Sol method.

        Parameters:
            struct_dir (str|Path):  Directory containing a finished VASP calculation on a structure.

            ignore_file (str):      Checks whether the provided file name exists in structure directory.
                                    Structure directories containing this file will return None.
                                    This parameter permits the filtration of structures that did not pass
                                    previous screening steps.
        
        Returns:
            Tuple[str, Dict]:       Tuple containing the name of the struct_dir and corresponding dict, 
                                    containing following useful data:
                                        - structure itself, 
                                        - its CHGCAR file (to modify charge density), 
                                        - its final energy (in eV/atom), used as E(N0).
        '''
        assert isinstance(struct_dir, PathLike)
        
        struct_dir: Path = Path(struct_dir)
        
        assert struct_dir.is_dir()

        files = set(file.name for file in struct_dir.iterdir())
        
        if ignore_file in files:
            return None

        struct_name  = struct_dir.name
        contcar_path = Path(struct_dir / 'CONTCAR')
        chgcar_path  = Path(struct_dir / 'CHGCAR')
        oszicar_path = Path(struct_dir / 'OSZICAR')
            
        structure    = Poscar.from_file(contcar_path).structure
        chgcar       = Chgcar.from_file(chgcar_path)
        final_energy_eV = Oszicar(oszicar_path).final_energy
        final_energy_eV_per_at = final_energy_eV / structure.num_sites

        struct_data = tuple(
            struct_name, 
            {
            'structure': structure, 
            'CHGCAR': chgcar, 
            'final_energy': final_energy_eV_per_at
            }
        )

        return struct_data

    def extract_all_vasp_data(
            base_dir: PathLike = '.', 
            ignore_file: str = 'rejected.txt', 
            workers: int = 1
    ) -> Dict[Dict[Structure, Chgcar, float]]:
        '''
        Extracts VASP data from a previous run for each structure directory in given directory.
        Keeps only data that are useful for Δ-Sol method.

        Parameters:
            base_dir (str|Path):    Directory containing structures subdirs to extract data from.

            ignore_file (str):      Checks whether the provided file name exists in each subdirectory.
                                    Structure directories containing this file will not be taken into account.
                                    This parameter permits the filtration of structures that did not pass
                                    previous steps.
        
        Returns:
            dict[dict]:             Dict with structure directory names as keys, 
                                    and a dict containing following data for corresponding structure as values:
                                        - structure itself, 
                                        - its CHGCAR file (to modify charge density), 
                                        - its final energy (in eV/atom), used as E(N0).
        '''

        assert isinstance(base_dir, PathLike)
        
        base_dir: Path = Path(base_dir)
        
        assert base_dir.is_dir()
        
        def is_directory(path: Path) -> bool:
            return path.is_dir()
        
        set_vasp_extractor = partial(extract_vasp_data, ignore_file=ignore_file)
        structs_dir_list   = list(filter(is_directory, base_dir.iterdir()))
        nbr_structs        = len(structs_dir_list)
        chunksize          = (min(nbr_structs // 100, 10) if nbr_structs >= 200 else 1)

        structs_data_list  = list(filter(
            process_map(
                set_vasp_extractor, 
                structs_dir_list, 
                workers=workers, 
                chunksize=chunksize, 
                desc='Extracting previous VASP outputs'
            )
        ))
        structs_data = dict(structs_data_list)

        return structs_data


    # Créer une fonction qui récupère les fichiers CHGCAR dans input_path, 
    # puis qui les convertit en objets CHGCAR dont on peut modifier la densité de charge
    # afin de les réécrire dans output_path avec les fichiers INCAR, POSCAR, POTCAR, KPOINTS.
    def chgcar_density_switch(chgcar: Chgcar, delta: float):
        '''
        Modifies provided Chgcar object's charge density by +/- delta, 
        then returns the two resulting Chgcar objects.
        '''
        e = 1.60217663e-19 # Coulomb.
        chgcar_plus, chgcar_minus = chgcar.copy(), chgcar.copy()
        chgcar_plus.data['total'] += delta*e
        chgcar_minus.data['total'] -= delta*e
        return chgcar_plus, chgcar_minus
    
    def struct_charge_switch(structure: Structure, delta: float):
        '''
        Modifies overall charge of provided structure by +/- delta, 
        then returns the two resulting structures.
        '''
        cation_struct, anion_struct = structure.copy(), structure.copy()
        cation_struct.set_charge(structure.charge + delta)
        anion_struct.set_charge(structure.charge - delta)
        return cation_struct, anion_struct

    # Réaliser les calculs VASP de E(N0 - n) et E(N0 + n).

    # Réutiliser la 1ère fonction pour récupérer les énergies des OSZICAR résultants.

    '''for struct in structures:
        N_0 = get_all_valence_electrons(struct)
        n   = get_delta_sol_el_ratio(struct, 'PBE', 'BEST')
        # E(N0) est à lire directement dans le fichier OSZICAR (dernier E0) de la relaxation.
        # Changer la valeur de la densité de charge dans CHGCAR :
        # Densité de charge ponctuelle n(r) = densité électronique ponctuelle rhô(r) x charge élémentaire e
        # Densité de charge CHGCAR nc(r) = densité de charge ponctuelle n(r) x Vol. maille V(maille)
        # Ccl: nc(r) = n(maille) représente la densité de charge par maille, on peut donc enlever/ajouter
        # une densité de charge n x e directement à nc(r) pour faire les calculs 
        # de minimization électronique, ioniquement statiques E(N0 - n) et E(N0 + n).
        # Les énergies seront également à récupérer dans le fichier OSZICAR.
        # Enfin, on pourra faire le calcul EFG = [E(N0 + n) + E(N0 - n) - 2E(N0)]/n.'''

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')


if __name__ == '__main__':
    main()