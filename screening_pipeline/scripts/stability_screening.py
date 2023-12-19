#!/usr/bin/python
'''
A script using previous VASP relaxations to compute the relative stability of given structures.
'''

########################################
# SYSTEM I/O MODULES

from typing import List, Tuple
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter
from pathlib import Path

########################################
# OPTIMIZATION MODULES

from itertools import count
from tqdm.contrib.concurrent import process_map

########################################
# PYTON MATERIALS GENOMICS PACKAGE

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import VaspInput

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:

    assert Path(args.input_dir).is_dir(), \
    f'{args.input_dir}: No directory found.'

    assert Path(args.output).is_dir(), \
    f'{args.output}: No directory found.'

    assert args.ignore.endswith('.txt'), \
    'Structure ignoring file must be a plain text file type (.txt).'

    assert args.workers >= 1, \
    'The number of workers cannot be negative or zero.'

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name          = 'stability_screening'
    prog_description   = 'A script using previous VASP relaxations to compute the relative \
                          stability of given structures.'
    prog_missing_steps = '''
        Missing steps to complete the script:
            - Determine the critical formation reaction of the structure
            - (Relax reference structures for calculations consistency)
            - Compare structure energy with the sum of reference structure energies
            - Aknowledge what to do with unstable structures
        '''
    helper_format = RawDescriptionHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_description, 
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
        '-i', 
        '--ignore', 
        type=str, 
        default='rejected.txt', 
        help='''Defines a file whose presence in a structure directory means it did not pass
        previous screening steps and should not be used in this calculation. 
        This file will also be written in structure directories that did not pass this step.
        WARNING: 
        If the name of this file is overwritten, care must be taken that it is the same file
        throughout every used screening steps to make sure rejected structures don't go further.''', 
        metavar='ignore_file.txt'
    )
    parser.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of parallel processes to spawn for parallelized steps.',
        metavar='int',
    )
    args: Namespace = parser.parse_args()

    assert_args(args)
    
    input_dir   = Path(args.input_dir)
    outdir      = Path(args.output)
    ignore_file = args.ignore
    workers     = args.workers


    # MAIN BLOCK

    from screening_pipeline.utils.vasp_io import batch_extract_vasp_data

    structs_data = batch_extract_vasp_data(
        method='convex_hull', 
        base_dir=input_dir, 
        ignore_file=ignore_file, 
        workers=workers
    )

    # Script préliminaire
    from pymatgen.core.composition import Element, Composition
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

    for name, data in structs_data:

        comp  = Composition(data['formula'])
        entry = PDEntry(
            composition=comp, 
            energy=data['final_energy'], 
            name=name
        )

        '''
        Fil directeur du premier script ci-dessous :

        0/ Définitions :

            1 - CH = Convex Hull (Diagramme de phases compositionnels avec énergies).

            2 - Grandeur d'une structure = nb d'éléments différents dans la structure.

        1/  Grouper les structures par CH indépendantes (Les structures ayant le plus d'éléments
            peuvent servir à définir les CH, ex : FeTiGe3O2 peut générer la CH Fe-Ti-Ge-O s'il 
            n'y a pas de structure plus grande contenant tout ces éléments).

            Pour faire ce groupage, il faut donc : 

            1 - Détecter les structures les plus grandes.

            2 - Si certaines ont au moins 2 éléments en commun, chercher les structures ayant une
                formule entièrement inclue dans le sous-enemble d'éléments communs.

            3 - S'il y en a, les grandes structures partiellement communes et les petites 
                correspondantes forment une unique CH. Sinon, les grandes structures formeront
                des CH séparées.
            
            4 - Ranger les structures dont la formule est inclue dans celle d'une des grandes 
                structures avec la grande structure correspondante.
            
            5 - S'il reste des structures non groupées, reprendre à l'étape 1 sur le sous-ensemble
                non-groupé.
        
        1.5/Récupérer les structures de références qui rentrent dans les groupes construits.

        2/  Fabriquer une CH par groupe de structures, en initialisant les éléments simples 
            à partir du contenu du groupe (donner l'attribut 'ref' aux structures de référence).
        
        3/  Calculer les ΔH dans chaque CH.

        4/  Eliminer les structures dont le ΔH est trop important.
        '''

        # Alternatives possibles :
        # - Chercher les matériaux du dump qui correspondent chimiquement.
        # - Construire le diagramme de phase correspondant à cet ensemble.
        # - Demander la décomposition et l'énergie /r à la CH du candidat.
        # - Ecrire le fichier de rejet si l'énergie dépasse le seuil autorisé.
        # Optimisation : rassembler d'abord tout les candidats rentrant dans une même CH, 
        # afin de ne construire chaque CH de référence utiles qu'une fois.
        # Ajouter les CH de références OQMD préconstruites dans un fichier utils.
        # Ecrire une fonction pour rassembler les structures selon la CH correspondante.
        # Ecrire une fonction qui compare en batch les structures d'un groupe avec sa CH, 
        # puis ajoute à leur data respectif la valeur 'delta_H'. 
        # Opti : Minimiser le nombre de CH à construire en groupant plus largement les candidats.
        # Opti : Matcher et éliminer les structures équivalentes à celles de référence 
        #        (StructureMatcher).
        # Opti : Fabriquer une CH dynamique à partir des candidats d'un groupe, 
        #        afin de voir les plus stables du groupe.
        # Opti : Construire les CH de référence expérimentales en réalisant des calculs statiques 
        #        sur les ~263k structures expérimentales de l'ICSD une seule fois (reproduction
        #        de l'article 00). On peut même élaguer en éliminant les structures exp contenant
        #        des gaz ou des terres rares (~184k restant).

    stop = datetime.now()
    print(f'elapsed time: {stop-start}')

if __name__ == '__main__':
    main()
