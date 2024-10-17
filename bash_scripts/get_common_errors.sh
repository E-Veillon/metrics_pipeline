#!/bin/bash
# Search for common slurm and VASP errors in slurm output and error files.

grep "error" $@
grep "Error" $@
grep "EEEEEEE" $@
grep "TIME LIMIT" $@
grep "99 F" $@
