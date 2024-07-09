#!/bin/bash
# The built-in bash "cp" command allows to copy several files and directories in one destination.
# This script is designed to do the reverse, i.e. distribute a single source file or directory to a range of destinations.
# If -r is used, the source directory and everything contained inside is distributed.
# Usage: ./multicopy.sh [-r] SOURCE DEST1 DEST2 ...

if [ "$1" == "-r" ]; then
	recursive=1
	shift
else
	recursive=0
fi

src_dir="$1"
shift

if [ $recursive -eq 1 ]; then
	for dest_dir in "$@"; do
		cp -r "$src_dir" "$dest_dir"
	done
else
	for dest_dir in "$@"; do
		cp "$src_file" "$dest_dir"
	done
fi

