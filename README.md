[Description général des objectifs de la pipeline]
...

[Définitions des métriques mesurées]
- Validity: No pair of atom in the structure are closer than 0.5 angstroms (50 pm).

- Stability, Uniqueness, Novelty (S.U.N.):

    - Stability: The structure's energy per atom above the convex hull of its chemical space is below a defined threshold (default 0.1 eV/atom).

    - Uniqueness: The structure is not equivalent to one previously encountered in the generation batch (Note: the first iteration of the structure is always considered unique, even if other structures are afterward considered equivalent to it).

    - Novelty: The structure is not equivalent to any structure used in the model training set.

- Average Root Mean Square Displacement (RMSD): Measure the mean squared distance between generated position and DFT equilibrium position of each ion in a structure, then compute the mean over all generated structures.

- Coverage (COV-P, COV-R):

    - Precision (COV-P):

    - Recall (COV-R):

- Fréchet ALIGNN Distance (FAD):

- Earth Mover's Distance (EMD) on density or energy:


[Fonctionnement détaillé des scripts]
0) Fichiers nécessaires avant de commencer :
    - Un fichier CIF contenant toutes les structures de référence utilisées comme input pour la génération avec votre modèle d'IA.
    - Un fichier CIF contenant toutes les structures générées à tester.

1) Preprocess.py - Ce script permet de trier et filtrer les structures générées selon plusieurs critères.
    Les critères de filtrage suivant sont activés par défaut, et élimineront du fichier de sortie les structures qui ne les respectent pas. Ils sont tous désactivables en passant différents drapeau dans la commande d'appel du script :

    - La structure est éliminée si elle contient des éléments chimiques de la famille des gaz rares ou des terres rares (désactivables respectivement avec les drapeaux --no-rare-gas-check et --no-rare-earth-check)

    - La structure est éliminée si elle n'est pas "valide" au sens de la métrique de Validité. Une structure est "valide" si elle ne contient aucune paire d'atome dont la distance est inférieure à 0.5 angstroms (désactivable avec le drapeau --no-valid-check)

    - La structure est éliminée si elle est équivalente à une structure précédente du fichier donné en entrée. Ce critère définit l'équivalence selon les paramètres de tolérance par défaut du StructureMatcher de pymatgen (désactivable avec le drapeau --no-equiv-match)

    Fonctions supplémentaires :

    - Les structures restantes sont triées selon leur composition chimique (non désactivable)

    - Par défaut, le script va chercher à obtenir le groupe de symétrie d'espace des structures pour l'afficher dans le fichier de sortie. Cette fonction ne sert qu'à la visualisation et n'est pas nécessaire au bon fonctionnement de la pipeline, donc si elle ne vous intéresse pas, désactivez la pour gagner en temps de calcul (désactivable avec le drapeau --no-symmetrization)

2) vasp_static_sun.py - Ce script fait appel au Vienna Ab-initio Simulation Package (VASP) afin de réaliser un calcul de minimisation électronique sur les structures sorties du filtrage de preprocess.py afin de mesurer leur énergie totale, nécessaire pour le calcul de stabilité de la métrique S.U.N. Les positions des ions dans les structures restent fixes lors de cette étape.

ATTENTION : Cette étape génère en sortie un répertoire par structure, chacun contenant plusieurs fichiers de sortie correspondant à un calcul VASP. Un grand nombre de répertoire peut être créé à l'emplacement spécifié pour la sortie (autant que de structures d'entrée), il est donc fortement recommandé d'allouer un répertoire vide uniquement dédié au stockage de ces données, qui peuvent également être relativement lourdes si certaines structures sont difficiles à calculer (i.e. nécessitent beaucoup d'itérations car éloignées de leur géométrie d'équilibre). De plus, une indexation doit être spécifiée pour traiter chaque structure indépendamment (l'index correspond à la position de la structure dans le fichier), aussi l'utilisation d'un tableau de jobs comme proposé dans le gestionnaire de jobs Slurm est également recommandé. Le script cif_counter.py permet de savoir efficacement le nombre exact de structures dans un fichier CIF, il peut être utile pour vous aider dans l'indexation.

3) phase_diag_energies.py - Ce script récupère les données générées par VASP via le script vasp_static_sun.py et les comparent à un fichier de base de données de structures connues afin de déterminer leur stabilité relative par rapport aux données connues du même système chimique via la construction d'enveloppes convexes, gérées par pymatgen. Ce script génère un fichier JSON (par défaut nommé "summary.json) contenant les informations de stabilité des structures générées (nom et chemin du répertoire de la structure, énergie au-dessus de l'enveloppe mesurée, la structure est-elle considérée stable ou pas).

4) vasp_relax.py - Ce script refait appel à VASP pour cette fois optimiser la géométrie des structures en entrée afin de trouver ses positions d'équilibre. Cette étape servira à mesurer la métrique de RMSD entre les structure générées et leur version calculée avec la Density Functional Theorie (DFT).

ATTENTION : Cette étape peut être très coûteuse en temps de calcul si les structures sont loin de leur position d'équilibre (jusqu'à plusieurs jours par structure avec 16 coeurs CPU). Il est recommandé de la réaliser sur un supercalculateur et de définir un walltime au bout duquel le calcul s'arrête et la structure est considérée comme trop loin de son équilibre pour être mesurée dans le RMSD. Certaines structures peuvent ne jamais trouver de position d'équilibre en finissant sur un erreur VASP, en atteignant le maximum d'itérations ioniques (optimisation des positions des ions) ou le walltime défini.

5) metrics.py - Ce script permet de calculer en une fois toutes les métriques supportées par la pipeline.
Une fois que tout les fichiers nécessaires ont été rassemblés ou générés par les scripts précédents, renseignez-les dans les arguments en ligne de commande correspondants. Si certaines métriques ne vous intéressent pas, vous pouvez désactiver leur calcul avec des drapeaux. Dans ce cas, les fichiers associés peuvent ne plus être nécéssaires et ne pas être chargés. Le script génère un fichier json listant les valeurs mesurées pour les différentes métriques, et mettra la valeur "null" sur les métriques désactivées.

ATTENTION :
- Pour mesurer la Coverage (Precision, Recall), vous avez besoin d'installer CrystalNN dans les dépendances.
- Pour mesurer la Fréchet ALIGNN Distance, assurez-vous d'avoir un accès internet là où le calcul des métriques est réalisé ou alors téléchargez le répertoire zippé du modèle ALIGNN préentrainé et déplacez-le à l'emplacement prévu par le programme.

[Workflow conseillé]