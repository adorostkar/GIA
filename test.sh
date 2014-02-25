#!/bin/bash
FILE=outc1i.txt
C=1

# Write ARGS directive
echo ARGS Refinements  Inner Threads  | tee -a ${FILE}
# Write caption depending on block size
if [ $C = "1" ]
then
    echo CAPTION:2-by-2 block structure | tee -a ${FILE}
else
    echo CAPTION:3-by-3 block structure | tee -a ${FILE}
fi

run experiments
for refine in 3 4 5 6 7
do
    for ompt in 1 4 8 16
    do
        for inner in 1e-1 1e-4 1e-8
        do
            echo Arg Refinements $refine Inner $inner Threads $ompt | tee -a ${FILE}
            export OMP_NUM_THREADS=$ompt
            ./elastic -r $refine -i $inner -c $C | tee -a ${FILE}
        done    
    done
done
