#!/bin/bash

##### DOCUMENTATION #####

SCRIPT_NAME="$(basename -- "${BASH_SOURCE[0]}")"

Help() {
cat << EOF

Usage: $SCRIPT_NAME <error_file1> [<error_file2> ...]

Search for common slurm and VASP errors in slurm output and error files.

Options
	-h, --help	Show this message and exit.

EOF
}

for arg in "$@"; do
	if [ $arg = "-h" -o $arg = "--help" ]; then
		Help
		exit 0
	fi
done

##### ACTUAL CONTENT #####

grep "error" $@
grep "Error" $@
grep "EEEEEEE" $@
grep "TIME LIMIT" $@
grep "99 F" $@
