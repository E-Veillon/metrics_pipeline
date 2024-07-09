#!/bin/bash

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

