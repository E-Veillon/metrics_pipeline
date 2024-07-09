#!/bin/bash

src_file="$1"
shift

for file in "$@"; do
	echo "diff $src_file and $file:"
	diff "$src_file" "$file"
done
