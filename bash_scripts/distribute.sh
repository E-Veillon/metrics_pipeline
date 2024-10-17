#!/bin/bash
# The built-in bash "cp" command allows to copy several files and directories in one destination.
# This script is designed to do the reverse, i.e. distribute a single source file or directory to a range of destinations.
# If -r is used, the source directory and everything contained inside is distributed.
# Usage: ./multicopy.sh [-r] SOURCE DEST1 DEST2 ...

declare -i src_defined=0
declare -i recursive=0
dest_dirs=()

while [ $# -gt 0 ]; do
    case "$1" in
        "-r"|"--recursive" )
	    recursive=1
	    echo "'recursive' flag detected"
	    shift 1
	    ;;
	-* )
	    echo "unrecognized argument specifier: $1"
	    exit 1
	    ;;
	* )
	    if [ $src_defined -eq 0 ]; then
		src="$1"
		echo "source file defined as $1"
		src_defined=1
	    else
		dest_dirs+=("$1")
		echo "added $1 to destinations"
	    fi
	    shift 1
    esac
done
echo "registered destinations: ${dest_dirs[@]}"
if [ $recursive -eq 1 ]; then
	for dest in ${dest_dirs[@]}; do
		cp -r "$src" "$dest"
	done
else
	for dest in ${dest_dirs[@]}; do
		cp "$src" "$dest"
	done
fi

