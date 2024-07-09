#!/bin/bash
# The built-in bash command "diff" allows to spot differences between 2 files.
# This script is designed to be a recursive approach of it, comparing a single source file with all destination files.
# Usage: ./multidiff.sh SOURCE DEST1 DEST2 ...
# Result:
# "diff SOURCE and DEST1:"
# {diff between SOURCE and DEST1}
# "diff SOURCE and DEST2:"
# {diff between SOURCE and DEST2} etc...

src_file="$1"
shift

for file in "$@"; do
	echo "diff $src_file and $file:"
	diff "$src_file" "$file"
done
