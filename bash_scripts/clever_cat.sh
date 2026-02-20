#!/bin/bash

##### DOCUMENTATION #####

SCRIPT_NAME="$(basename -- "${BASH_SOURCE[0]}")"

Help() {
cat << EOF

Usage: $SCRIPT_NAME [--sep SEPARATOR] [--target NEWFILE] <file1> <file2> [<file3> ...]

More convenient 'cat' command to concatenate data from text files.

Options
	-h, --help		Show this message and exit.

	--sep SEPARATOR		Give alternative string that should be used to separate each file's data.
				Defaults to a single newline character.

	--target NEWFILE	Name of the concatenated file. Defaults to "concatenated.txt".

EOF
}

##### ARGUMENTS PROCESSING #####

declare -i isFlag=0
declare -i isSep=0
declare -i isTarget=0
declare -a files=()
sep=""
target="concatenated.txt"
Counter=0
firstFile=""

CheckFlag () {
	isFlag=0
	if [ $1 = "-h" -o $1 = "--help" ]; then
		Help
		exit 0
	elif [ $1 = "-s" -o $1 = "--sep" ]; then
		isFlag=1
		isSep=1
	elif [ $1 = "-t" -o $1 = "--target" ]; then
		isFlag=1
		isTarget=1
	fi
}

for arg in "$@"; do
	CheckFlag "$arg"
	if [ $isFlag -eq 1 ]; then
		continue
	elif [ $isSep -eq 1 ]; then
		sep="$arg"
		isStep=0
	elif [ $isTarget -eq 1 ]; then
		target="$arg"
		isTarget=0
	else
		files+=("$arg")
	fi
done

##### ACTUAL CONTENT #####

mkdir -p ".cat_tmp"

for file in ${files[@]}; do
	if [ $Counter -eq 0 ]; then
		firstFile="$file"
	else
		tmp=".cat_tmp/cat_tmp_${Counter}.txt"
		echo "$sep" | cat "$firstFile" - "$file" > "$tmp"
		firstFile="$tmp"
	fi
	Counter=$(($Counter + 1))
done

mv "$tmp" "$target"
rm -r ".cat_tmp"
