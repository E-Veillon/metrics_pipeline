# Poetry Cheat sheet


Initialiser un nouveau paquet avec poetry :
```bash
poetry new [NomDuPaquet]
```

Lancer un script dans l'environnement virtuel (sans extension) :
```bash
poetry run [NomDuScript]
```

Mettre à jour les dépendances et les scripts pris en compte :
```bash
poetry install
```

Mettre à jour le fichier poetry.lock :
```bash
poetry lock
```

Utilisation de la pipeline sur Jean Zay :

1) Création de l'environnement de calcul :
    1.1 - Créer dans $SCRATCH (meilleure partition pour la lecture-écriture) un répertoire de stockage des données 
          pour un batch de structures générées.
    1.2 - Mettre dans le répertoire le fichier CIF concatené de la génération.

2) pré-processing :
    2.1 - Dans le script de soumission "preproc_job.slurm", vérifier et modifier si besoin :
        * Le nom du répertoire du batch dans la ligne "export RUNDIR="
        * Le nom du fichier CIF concatené dans la ligne "export INFILE="
        * Le nom voulu pour le fichier CIF de sortie dans la ligne "export OUTFILE="
        * Les flags à utiliser dans la ligne d'exécution ("poetry run symmetrize --help" pour la liste des flags)

    2.2 - Soumettre le script preproc_job.slurm -> le fichier CIF avec le nom voulu est créé à côté du premier.

3) Relaxation :
    3.1 - Dans le script de soumission "relax_job.slurm", vérifier et modifier si besoin :
        * Le nom du répertoire du batch dans la ligne "export RUNDIR="
        * Le nom du fichier d'entrée qui correspond à la sortie du pré-processing dans la ligne "export INFILE="
    
    3.2 - Utilisation des Job Arrays :
        * Vérifier dans la ligne "#SBATCH --array=" du script de soumission que le nombre de jobs correspond 
          au nombre de structures dans le fichier d'entrée.
        * Le script cif_counter.py dans screening_pipeline/scripts peut être utilisé pour compter rapidement 
          les structures dans le CIF avec la commande "python cif_counter.py chemin/du/fichier.cif"

        ATTENTION : Par défaut, les Job Arrays ne peuvent accepter que 1000 jobs d'un seul coup. S'il y a plus de 1000 structures à relaxer, 
        il faudra soumettre plusieurs jobs avec des tranches d'indices différents, car le job N ira chercher la Nième structure du fichier.
        (exemple : pour 1500 structures, il faudra d'abord soumettre un array allant de 0 à 999, puis un second array allant de 1000 à 1499)

    3.3 - Soumettre le script relax_job.slurm -> un répertoire "Relaxations" est créé dans le répertoire du batch, 
          dans lequel chaque structure possède un sous-répertoire de calcul VASP.

4) Stabilité :
    4.1 - Dans le script de soumission "stability_job.slurm", vérifier et modifier si besoin :
        * Le nom du répertoire du batch dans la ligne "export RUNDIR="

    4.2 - Soumettre le script stability_job.slurm -> un fichier "summary.json" est créé dans le répertoire "Relaxations", qui donne l'énergie
          par rapport à l'enveloppe convexe des structures issues du dataset de chaque structure générée, et indique si elle est considérée
          stable ou non.
