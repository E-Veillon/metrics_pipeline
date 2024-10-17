#!/bin/bash
# The built-in bash command "diff" allows to spot differences between 2 files.
# This script is designed to be a recursive approach of it, comparing a single source file with all destination files.
# Equivalent to "diff --from-file=SOURCE DEST1 DEST2 ...", but with a title line naming compared files before each comparison.
# Usage: ./multidiff.sh SOURCE DEST1 DEST2 ...

src_file="$1"
shift

for file in "$@"; do
	echo "diff $src_file and $file:"
	diff "$src_file" "$file"
done
