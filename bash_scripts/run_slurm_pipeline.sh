#!/bin/bash

##### DOCUMENTATION #####

SCRIPT_NAME="$(basename -- "${BASH_SOURCE[0]}")"

Help() {
cat << EOF
Usage: $SCRIPT_NAME <script1.slurm> [<script2.slurm> ...]

Submit slurm jobs as a linear pipeline where each job script is run whenever the previous one succeeded.
Scripts will be submitted in the order given by passed arguments.
A failed run for one script will end the pipeline.

Options
	-h, --help	Show this message and exit.

EOF
}

for arg in "$@"; do
	if [ $1 = "-h" -o $1 = "--help" ]; then
		Help
		exit 0
	fi
done

##### ACTUAL CONTENT #####

for script in "$@"; do
	if [ -z $jobid ]; then
		jobid=$(sbatch "$script" | cut -d " " -f 4)
	else
		jobid=$(sbatch --dependency=afterok:"$jobid" "$script" | cut -d " " -f 4)
	fi
	echo "Submitted batch Job $jobid for script $script"
done

