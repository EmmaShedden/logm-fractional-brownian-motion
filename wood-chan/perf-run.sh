#!/bin/bash
# q-bound : int
# num-repeats : int
# output-filename : string

single_test () {
	# $1: q
	# $2: output filename
	echo -n "$1, " >> $2
	perf stat -x , python3 main.py 2>> $2
	echo "$1 complete"
}

if [ $# -eq 3 ]
then
	# run each q value n times
	for (( i=0;  i<$1; i++ )); do
		for (( j=0; j<$2; j++ )); do
			#echo "$1 $2 $3 $i $j";
			single_test $i $3;
		done
	done
else
	echo "missing arguments: ./perf-run.sh <q-bound> <num-repeats> <output-filename>"
fi

