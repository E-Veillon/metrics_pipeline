#!/bin/bash

WarningCounter () {
	warn_count=$(grep -c "$2" $3)
	if [ $warn_count -gt 0 ]; then
		echo "- $1 warnings: $warn_count"
	fi
}

FinishedOnError () {
	has_error="no"
	error_count=$(grep -c "EEEEEEE" "$1")
	if [ $error_count -gt 0 ]; then
		has_error="yes"
	fi
}

declare -i cumul_warn_count=0

echo " "
echo "===== BEGINNING OUTCAR PARSING ====="

for file in "$@"; do
	echo " "
	echo "------------------------------"
	echo " "
	echo "$file:"
	echo "Last iteration:"
	grep Iteration $file > "grep_iter.tmp"
	tail -1 "grep_iter.tmp" > "tail_iter.tmp"
	echo "Ionic: $(cut -c 53-56 "tail_iter.tmp")"
	echo "Electronic: $(cut -c 58-61 "tail_iter.tmp")"
	echo "Last total energy:"
	grep TOTEN $file > "grep_toten.tmp"
	tail -1 "grep_toten.tmp"
	grep "energy(sigma->0)" $file > "grep_sigma_0.tmp"
	tail -1 "grep_sigma_0.tmp"
	total_warn_count=$(grep -c "WW  WW" $file)
	echo "Number of Warnings: $total_warn_count"
	if [ $total_warn_count -gt 0 ]; then
		WarningCounter "Interatomic distances" "distance between some ions" $file
		cumul_warn_count=$(($cumul_warn_count + $warn_count))
		WarningCounter "Electronic not converged" "electronic self-consistency" $file
		cumul_warn_count=$(($cumul_warn_count + $warn_count))
		WarningCounter "Lattice symmetry class" "reciprocal lattice and k-lattice" $file
		cumul_warn_count=$(($cumul_warn_count + $warn_count))
		WarningCounter "Fermi occu variations" "Tetrahedron method does not include" $file
		cumul_warn_count=$(($cumul_warn_count + $warn_count))
		WarningCounter "Lattice len over 50" "lattice vectors is very long (>50 A)" $file
		cumul_warn_count=$(($cumul_warn_count + $warn_count))
		warns_left=$(($total_warn_count - $cumul_warn_count))
		if [ $warns_left -gt 0 ]; then
			echo "- other warnings: $warns_left"
		fi
	fi
	FinishedOnError $file
	echo "Finished on VASP Error: $has_error"
	if [ "$has_error" == "yes" ]; then
		grep -B 2 -A 20 -m 1 "EEEEEEE" $file
	fi
done

rm "grep_iter.tmp"
rm "tail_iter.tmp"
rm "grep_toten.tmp"
rm "grep_sigma_0.tmp"

echo " "
echo "===== PARSING ENDED SUCCESSFULLY ====="
echo " "
