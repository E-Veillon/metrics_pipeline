metrics_pipeline v1.2.3

The main goal of this pipeline is to measure some of the metrics that are often used in recent works that propose generative AI models to predict new materials from an input of known materials.

Whatever the processing done to the data before and after passing it to the model, the input and generated materials that will be used in the pipeline for evaluation have to be described, at least, by unit cell parameters, atomic compositions and positions inside the unit cell, in a CIF file format to make the pipeline work properly.

# Installation

We recommend using poetry to install dependencies.
Install poetry from PyPI:
```bash
pip install poetry
```
Make it create a virtual environment using provided pyproject.toml once positionned in this directory:
```bash
poetry install
```

# Usual workflow to evaluate metrics

## 0 - Needed data before beginning

- A CIF file containing all reference structures used as input to the model.

- A CIF file containing all generated structures to test.

- A JSON file containing following minimal informations about all known structures in considered database or combination of databases, used to build reference convex hulls in order to measure stability: a unique ID for each data in a "entry_id" key, a full formula of all atoms in the unit cell (e.g. a carbon diamond cell must be given as C8 and not only C) in a "composition" key, and a total energy in electronvolts (eV) in a "final_energy" key.

WARNING: Note that material databases often give total energies in electronvolts per atom (eV/atom), but when defining the data points, the program assumes total energies are in eV and divides them by the number of atoms in the unit cell (hence the need for a non-reduced formula), thus to have the correct data, it can be necessary to modify manually downloaded data to correspond to what is needed.

NOTE: The "dataset_loader.py" script can automatically download and format entries from the Materials Project using their API, and format an already downloaded JSON file of Open Quantum Materials Database (OQMD) entries if the right fields are included (see script documentation with the "python dataset_loader.py --help" command for more details). For other databases, the OQMD JSON formatting can still work if the fields have the same name and contain the same information as OQMD fields.

## 1 - Preprocess.py

Generate 3 output files from the same "generated structures" file input:

- one containing valid structures with all other features deactivated, for Validity metric;

- one containing unique structures by deactivating all features other than the comparator, for S.U.N. metric's Unicity part (see below for metrics definitions);

- one containing valid AND unique structures only for further processing.
    
Note: if rare elements must be discarded from the generated data, for consistency they must be also discarded from all reference files to not get wrong percentage results on metrics measurements (for the model input CIF file) and save time on filtering unused structures (for the JSON database file).

## 2 - vasp_static_sun.py

This script calls on the Vienna Ab-initio Simulation Package (VASP) to compute total energies of structures with DFT. First, prepare an empty directory for output data as this step generates a lot (one directory per structure computed, containing several configuration and output files each). Pass the CIF file containing only valid and unique structures to it, the path to the prepared output directory to the "--output" argument and an index to the "--task-index" argument to tell which structure in the file must be calculated.

NOTE: This step needs parallelisation to not take a prohibitive computation time. It is higly recommended to run it on a supercalculator. In this case, taking 16 CPU cores per structure works fine in general. If a job management software that proposes job arrays is used (like slurm), it is recommended to use it and pass the sub-job array ID in the "--task-index" argument to treat all (or a big part of) the structures in one go. 

## 3 - phase_diag_energies.py

Once ALL calculations from the previous step are finished in the output directory, pass the path to it to this script with the database JSON file in the "--reference" argument to compute reference convex hulls that corresponds to the generated structures chemical spaces and compare them to the hulls to determine their phase stability for S.U.N. metric "Stability" part (see below for metrics definitions).
A single JSON summary file will be generated at the end, containing for each structure its directory name and path, its energy per atom above the hull, and whether it is considered stable according to the threshold set in the "--limit" argument (0.1 eV/atom by default).

## 4 - vasp_relax.py

This script is only used to compute Average RMSD metric and is the most time consuming one.
The VASP software is called to optimise given structures geometry.
As for vasp_static_sun.py, it is prohibitively time consuming if not parallelised. See step 2) for parallelisation recommendations on a supercalculator.
The CIF file containing valid and unique structures generated from step 1 or the JSON file generated from step 3 can be passed. As only structures that did not issued VASP errors and converged normally in the "simple" static calculation are present in the JSON summary file, it is recommended to pass this one in order to limit time consuming errors during this step. However, this method uses the paths indicated in the JSON file to collect structures data, hence make sure when running this step that structures data are still in the right location, and if not it is possible to update the location stated in the file by using the "summary_paths_update.py" script. This script modifies all paths in the JSON file to a new path, hence make sure that all structure data stated in the JSON file are still present in the same directory to avoid "FileNotFound" errors.

## 5 - metrics.py

Finally, once all data have been collected, pass all necessary files in corresponding arguments to this script to compute the metrics. If not all implemented metrics interest you, flags can be passed to deactivate the computation of each metric individually. Files that are no longer needed for deactivated metrics do not have to be given and are skipped if given anyway. A single JSON file will be generated, containing all metrics values, with a "null" for deactivated ones.

NOTE:
- To measure Coverage metrics (Precision, Recall), you need to have the CrystalNN model python software installed in dependencies.
- To measure Fréchet ALIGNN Distance (FAD), make sure that the processor units computing this script have an internet access (often not the case for supercalculator nodes), or download beforehand the zipped directory containing the pretrained ALIGNN model at https://figshare.com/ndownloader/files/31458811 and place it at the location "~/.cache/materials-toolkit/pretrained/alignn/mp/e_form.zip" (at the moment only Materials Project training for ALIGNN was used and tested with this pipeline).

----------------------------------------

# Definitions of measured metrics

- **Validity:** No pair of atom in the structure are closer than 0.5 angstroms (50 pm).
The metric value is the percentage of valid structures inside given batch.

- **Stability, Uniqueness, Novelty (S.U.N.):**

    - *Stability:* The structure's energy per atom above the convex hull of its chemical space is below a defined threshold (default 0.1 eV/atom).

    - *Uniqueness:* The structure is not equivalent to one previously encountered in the generation batch (Note: the first iteration of the structure is always considered unique, even if other structures are afterward considered equivalent to it).

    - *Novelty:* The structure is not equivalent to any structure used in the model training set.

The final metric value is the percentage of structures validating the 3 conditions inside given batch.
For information, several percentages are also measured for sub-combinations of the conditions.

- **Average Root Mean Square Displacement (RMSD):** Measures the mean squared distance between generated position and DFT equilibrium position of each ion in a structure, then computes the mean over all structures in given batch.

- **Coverage (COV-P, COV-R):**
 
    - *Precision (COV-P):* Proportion the generated materials getting in the ground truth materials distribution. In other words, how many generated materials are of high quality.

    - *Recall (COV-R):* Proportion of the ground truth materials getting in the generated materials distribution. In other words, how many ground truth materials are correctly predicted.

More details about Coverage metrics in the work of Ganea et al:
Octavian-Eugen Ganea, Lagnajit Pattanaik, Connor W Coley, Regina Barzilay, Klavs F Jensen,
William H Green, and Tommi S Jaakkola. Geomol: Torsional geometric generation of molecular
3d conformer ensembles. arXiv preprint arXiv:2106.07802, 2021. 8, 18

- **Fréchet ALIGNN Distance (FAD):**
Compare probability distributions of ground truth and generated materials using ALIGNN model's predictions.

- **Earth Mover's Distance (EMD) on density or energy:**
Compare CrystalNN generated fingerprints distributions of ground truth and generated materials against a given property.

----------------------------------------

# Supplementary optional scripts

- *Band Gap evaluation*
    It is possible to estimate the fundamental band gaps of the materials with "vasp_static_dsol.py" and "dsol_bandgap.py" scripts. It uses the Δ-Sol method proposed by Chan et al.
    This method performs at least 3 electronic minimizations with different amount of electrons considered in the material, then computes a ponderated difference in total energy between them to estimate the gap.
    More infos about the method in the original paper referenced below.
    If uncertainty estimation is needed, the "--with-uncertainties" flag can be given, then 7 computations per structure are performed instead of 3.

    *WARNING:* Uncertainty calculations reliability was not rigorously tested, the feature may sometimes lead to weird results. Future updates on the pipeline may try to correct this behaviour.

    *WARNING 2:* This feature uses VASP to calculate band gaps, if you need to use it, it is highly recommended to first prepare an empty output directory to allow VASP to write data inside (one directory for each evaluated structure will be written, itself containing directories for each VASP calculation on the structure).

    *WARNING 3:* Δ-Sol method does not support structures containing f-block elements.

    Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    
- *Structure perturbation*
    This feature allows efficiently generating randomly perturbed materials from ones given as input.
    In this random generation, atomic positions and each cell parameter can be perturbed in a user given range.
    If a generated material gets a too weird shape or volume that makes it obviously unphysical, it issues a warning, discards the unphysical data and retry a random generation until it respects minimal physicality conditions or there is no retry attempts left.

- **Bash scripts**
    Some convenient little bash scripts to help data management between calculation steps.

    - The *"distribute.sh"* script allows to distribute a single file or directory (through the use of the -r flag) to several locations at once. It can be particularly convenient to distribute a same input configuration to several VASP calculations efficiently without risk of errors for example.

    - The *"get_common_errors.sh"* script can be used on several slurm output and error files (for VASP calculations steps parallelized on a supercalculator using slurm as job management software particularly) to catch slurm and VASP errors that were eventually issued in their most common form. This script is basically just a chain of "grep" bash commands applied on common error patterns, it is not guaranteed to catch all errors if they happen to not follow the usual patterns.

    - The *"multidiff.sh"* script is a simple extension of the bash "diff" command, that basically does the same as the "diff --from-file=input_file output_files" bash command line, but with a header at the beginning of each comparison clearly stating which files are compared.

    - The *"outcar_parser.sh"* script allows a quick and compact parsing of several OUTCAR files' main informations (Total number of ionic and electronic steps done, last printed total energies, issued VASP warnings and errors). By default the parsing is printed on terminal, but if there are a lot of OUTCAR files to parse at once (> 10), it is recommended to output the command in a file. if this script is used on anything else than a OUTCAR file, there is no error check to verify it, therefore it may make a weird output for this file.

    - The *"spg_mapping.sh"* script allows a comprehensive parsing of spacegroup numbers distributions in symmetrized concatenated CIF files. The *"_symmetry_Int_Tables_number"* and *"_space_group_IT_number"* line headers are searched to define spacegroup number of all structures in the screened file. The former one is automatically used by Pymatgen when writing CIF files, such as the ones that are outputed by the *preprocess.py* script of the pipeline.

An emphasis is made on the fact that all scripts described in this section are purely optional, and do not affect the data or quality of the main metrics evaluation pipeline in any way. They are not guaranteed to work properly under unusual conditions or in some specific calculations cases, therefore if they are used, their outputs must be considered with care. Future updates may test more thoroughly these features and correct them if necessary.
