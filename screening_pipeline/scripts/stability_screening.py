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
    from typing import Sequence, Any, Union, Iterable
    from pymatgen.core.composition import Composition, Element
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    #from itertools import product, combinations
    from screening_pipeline.utils.data_process import group_by_dim

    '''def group_by_convex_hull(structs_data: dict) -> List[List[dict]]:
        
        'Group every structure in groups to make convex hulls from.'
        
        assert isinstance(structs_data, dict) and len(structs_data) > 0, \
        f'Invalid input provided ({structs_data}), it either was not a dict or was an empty dict.'

        groups    = []
        comp_list = [(name, data['composition'], len(data['composition'])) for name, data in structs_data]
        
        while comp_list:

            max_size          = max([tup[2] for tup in comp_list])
            max_sized_structs = list(filter(lambda t: t[2] == max_size, comp_list))
            smaller_structs   = list(filter(lambda t: t[2] < max_size, comp_list))

            if len(max_sized_structs) == 1:
                # All structures whose composition is only made of elements 
                # that are in the bigger composition are added in the group 
                # and removed from the main list.
                max_sized_comp = max_sized_structs[0][1]

                included_small_comps = list(filter(
                    lambda t: all(elt in max_sized_comp for elt in t[1]), smaller_structs
                ))
                group = [*max_sized_structs, *included_small_comps]
                groups.append(group)
                for struct in group: comp_list.remove(struct)
                continue

            pairwise_compare = list(combinations(max_sized_structs, 2))

            for comp1, comp2 in pairwise_compare:
                # For each possible pair of big struct, elements in common are searched
                common_elts    = list(filter(lambda elt: elt in comp2, comp1.elements))
                nb_common_elts = len(common_elts)

                if nb_common_elts < 2: continue

                # Vérifier la composition commune avec les autres compositions mères
                # plusieurs cas :
                #   - le pattern ne correspond à aucune autre structure
                #   - le pattern correspond à une seule autre structure
                #   - le pattern correspond à plusieurs autres structures

                # All structures whose composition is only made of elements 
                # that are in the common composition are added in the group 
                # and removed from the main list, along with the 2 bigger ones.
                common_comp          = Composition(common_elts)
                candidates           = list(filter(lambda t: t[2] <= nb_common_elts, smaller_structs))
                included_small_comps = list(filter(lambda t: all(elt in common_comp for elt in t[1]), candidates))
                group = [comp1, comp2, *included_small_comps]
                groups.append(group)
                # Manque des choses ici...'''

    def calculate_instability_energies(structs_data: dict):
        '''
        Construct an adaptive convex hull for each structure according to their composition.
        A binary structure does not need comparison with higher order structures.
        For a higher order structure, smaller convex hulls can be combined into one of correct 
        composition to put the structure in.

        Parameters:
            structs_data (dict):    A dict containing following data about each structure:
                                        - structure directory name (dict's keys), 
                                        - the Structure object, 
                                        - its composition, 
                                        - its relaxed energy (in eV).
        Returns:
            Dict: The same data with all ΔH calculated in 'delta_H' keys.
        '''
        
        assert isinstance(structs_data, dict) and len(structs_data) > 0, \
        f'''Invalid input provided, it either was not a dict or was empty.
            Detected type: {type(structs_data)}.
            Detected length: {len(structs_data)}.'''

        groups = group_by_dim(structs_data)
        # A ce stade, groups = [
        #                       [], 
        #                       [Elements], 
        #                       [[Binary 1 (ex Fe-O)], [Binary 2 (ex Mn-O)], ...], 
        #                       [[Ternary 1 (ex Fe-Mn-O)], [Ternary 2 (ex Fe-Co-O)], ...], 
        #                       ...
        #                      ]
        # On peut donc commencer à construire les Convex Hulls.
        cached_pds: dict[str, dict] = {}

        for elts_nbr in range(2, max_elts_nbr + 1):
            for comp_region in groups[elts_nbr]:

                elements   = list(filter(lambda elt: elt in comp_region[0].elements, groups[1]))
                pd_name    = '-'.join([elt.symbol for elt in elements])

                #TODO: finish cached diagrams conditional block
                if elts_nbr > 2 and 'previous diagrams are included in this one':

                    comp_pd = init_pd_from_cache(
                        ref_elts=pd_name, 
                        cached_pd_data=cached_pds, 
                        new_data=comp_region
                    )

                else:

                    comp_pd = init_pd_from_scratch(
                        ref_elts=pd_name, 
                        structs_data=comp_region
                    )

                for entry in comp_pd.all_entries[elts_nbr:]:
                    delta_H = comp_pd.get_e_above_hull(entry)
                    structs_data[entry.name]['delta_H'] = delta_H

                cached_pds[pd_name] = comp_pd.as_dict()
        
        return structs_data
    
    def get_elements(
            elts_data: str|Iterable[str|int|Element]
        ) -> List[Element]:
        '''
        Flexible converter to get a list of unique Element objects from a single string or any 
        iterable providing valid element symbols or atomic numbers, or a mixture of the two.

        Parameters:
            elts_data (str|[str|int|Element]):  The data to parse Elements objects from.

                                                If a single string is provided, it can either 
                                                be a raw formula (eg. 'FePO4') or a composition 
                                                string containing element symbols separated by 
                                                '-' (eg. 'Fe-P-O').

                                                If an iterable is given, it can contain valid 
                                                element symbols, atomic numbers and/or Element 
                                                objects.

        Raises: 
            ValueError if some of the given data does not represents valid elements.

        Returns: 
            A list of parsed Element objects.
        '''

        assert isinstance(elts_data, Iterable)

        if isinstance(elts_data, str):
            elts_list = Composition(''.join(elts_data.split(sep='-')), strict=True).elements
        
        else:
            assert all([isinstance(elt, (str, int, Element)) for elt in elts_data])
            elts_list = Composition([(elt, 1) for elt in elts_data], strict=True).elements

        return elts_list

    def init_pd_from_cache(
            ref_elts: Union[str, Iterable], 
            cached_pd_data: dict, 
            new_data: dict[str, dict]|Sequence[Tuple]
        ) -> PhaseDiagram:
        '''
        Initialize a higher order PhaseDiagram object by combining data from lesser ones already computed.
        Avoid expensive construction of complex phase diagrams and redundance in computations.
        If some parts of the higher diagram are not computed yet, this function computes them

        Parameters:
            ref_elts (str|Iterable):        The  elemental references of the new phase diagram.
                                            
                                            If a single string is provided, it can either 
                                            be a raw formula (eg. 'FePO4') or a composition 
                                            string containing element symbols separated by 
                                            '-' (eg. 'Fe-P-O').

                                            If an iterable is given, it can contain valid 
                                            element symbols, atomic numbers and/or Element 
                                            objects.
            
            cached_pd_data (dict):          A dict referencing computed diagrams by their name, 
                                            and containing their MSONable dicts, as returned by
                                            PhaseDiagram.as_dict() method.
            
            new_data (dict|[Tuple]):        Structure data to put in the phase diagram in addition to 
                                            previous diagrams, typically structures of same elemental
                                            composition as the new diagram itself. 
                                            They must be provided as a dict of {name: data} or a sequence 
                                            of (name, data) tuples, where 'name' is the name of the 
                                            structure directory where data were took from, and 'data' is 
                                            a dict containing same infos as provided by the function
                                            'screening_pipeline.utils.vasp_io.extract_vasp_data_for_convex_hull'.
        
        Returns:
            The constructed PhaseDiagram object.
        '''

        assert isinstance(ref_elts, (str, Iterable))
        assert isinstance(cached_pd_data, dict)
        assert isinstance(new_data, (dict, Sequence))

        ref_elts    = get_elements(ref_elts)
        new_pd_dim  = len(ref_elts)
        new_pd_data = {
            "@module": PhaseDiagram.__module__, # OK
            "@class": PhaseDiagram.__name__, # OK
            "all_entries": [], # à calculer
            "elements": [elt.as_dict() for elt in ref_elts], # OK
            "computed_data": {
                "facets": None, # utiliser get_facets une fois tout les points dans le tableau
                "simplexes": None, # transformer les facets en Simplex et le lister ici
                "all_entries": [], # à calculer
                "qhull_data": None, # numpy.ndarray, à caster en liste et remettre en array ensuite
                "dim": new_pd_dim, # OK
                "el_refs": [ # OK
                    (elt, PDEntry(
                        composition=Composition(elt), 
                        energy=0.0, 
                        name=elt.symbol, 
                        attribute='element_ref'
                    )) for elt in ref_elts
                ], 
                "qhull_entries": [] # à calculer
            }
        }

        all_sub_pd_list = list(filter(
            lambda pd_name: all(elt in ref_elts for elt in get_elements(pd_name)), 
            cached_pd_data.keys()
        ))

        for sub_dim in reversed(range(2, new_pd_dim)):
            dim_sub_pd_list = list(filter(
                lambda sub_pd: len(get_elements(sub_pd)) == sub_dim, 
                all_sub_pd_list
            ))
            if not dim_sub_pd_list: continue
            # TODO: Il faut éviter d'ajouter des sous-diagrammes si les arêtes 
            #       correspondantes sont déjà satisfaites.
            new_pd_data['all_entries'] += list(set(sum(
                cached_pd_data[pd_name]['all_entries'] for pd_name in dim_sub_pd_list
            )))


        new_pd_data['computed_data']['all_entries'] = [
            PDEntry.from_dict(entry) for entry in new_pd_data['all_entries']
        ]

        
        '''new_pd_data['all_entries'] += [
            PDEntry(
                composition=Composition(data['structure']), 
                energy=data['final_energy'], 
                name=name, 
                attribute='generated'
            ).as_dict() for name, data in (new_data.items() if isinstance(new_data, dict) else new_data)
        ]

        new_pd_data['all_entries'] = list(set(new_pd_data['all_entries']))
        new_pd_data['elements']    = list(set(new_pd_data['elements']))
        
        new_pd_data['computed_data']['all_entries'] = [
            PDEntry.from_dict(entry) for entry in new_pd_data['all_entries']
        ]

        new_pd_data['computed_data']['el_refs'] = [
            (Element.from_dict(elt), PDEntry(
                composition=Composition(Element.from_dict(elt)), 
                energy=0.0, 
                name=Element.from_dict(elt).symbol, 
                attribute='element_ref'
            )) for elt in new_pd_data['elements']
        ]'''

        '''
        diagram1 = cached_pd_data['diagram1']
        diagram2 = cached_pd_data['diagram2']
        ...

        entry_list = []

        for struct in new_data:
            entry = PDEntry(
                Composition(struct), 
                energy=struct.final_energy, 
                name=struct.name, 
                attribute='generated'
            )
            entry_list.append(entry)

        computed_data = {
        "@module": PhaseDiagram.__module__, 
        "@class": PhaseDiagram.__name__, 
        "all_entries": diagram1["all_entries"] + diagram2["all_entries"] + ... + entry_list, 
        "elements": diagram1["elements"] + diagram2['elements'] + ..., 
        "computed_data": diagram1["computed_data"] + diagram2["computed_data"] + ...
        }
        '''
        new_pd = PhaseDiagram.from_dict(dct=new_pd_data)
        
        return new_pd

    def init_pd_from_scratch(
            ref_elts: Union[str, Iterable], 
            structs_data: dict[str, Any]|Sequence[Tuple[str, Any]]
        ) -> PhaseDiagram:
        '''
        Compute a PhaseDiagram object from given elements and structure data.

        Parameters:
            ref_elts (str|Iterable):        The  elemental references of the new phase diagram.
                                            
                                            If a single string is provided, it can either 
                                            be a raw formula (eg. 'FePO4') or a composition 
                                            string containing element symbols separated by 
                                            '-' (eg. 'Fe-P-O').

                                            If an iterable is given, it can contain valid 
                                            element symbols, atomic numbers and/or Element 
                                            objects.

            structs_data (dict|Sequence):   The structures data to put into the diagram.
        
        Returns:
            The constructed PhaseDiagram object.
        '''

        ref_elts = get_elements(ref_elts)

        if isinstance(structs_data, dict):
            structs_data = list(structs_data.items())

        entry_list = [
            PDEntry(
                composition=Composition(elt), 
                energy=0.0, 
                name=elt.symbol, 
                attribute='element_ref'
            ) for elt in ref_elts
        ] + [
            PDEntry(
                composition=struct[1]['composition'], 
                energy=struct[1]['final_energy'], 
                name=struct[0], 
                attribute='generated'
            ) for struct in structs_data
        ]

        new_pd = PhaseDiagram(
            entries=entry_list, 
            elements=ref_elts.elements
        )

        return new_pd

        # Les cristaux purs devraient être une énergie de référence.
        # Les binaires doivent se voir assigner une CH 1D.
        # Les ternaires doivent se voir assigner une combinaison de 3 CH 1D si possible, 
        # sinon 2 CH 1D + 1 ligne vide entre les deux élts non-reliés, 
        # sinon 1 CH + 1 elt relié aux 2 références par des lignes vides, 
        # sinon une CH 2D vierge avec ses 3 elts pour références.
        # Même principe pour les structures d'ordre supérieur.

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
