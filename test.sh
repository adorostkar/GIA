#!/bin/bash

for refine in 2 3 4 5 6
do
	for inv in 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8
	do
		for schur in 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8
		do
			echo ./lib/elastic-2d -r $refine -i $inv -s $schur | tee -a out.txt
			./lib/elastic-2d -r $refine -i $inv -s $schur | tee -a output.out
		done
	done
done