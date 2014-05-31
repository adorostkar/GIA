#!/bin/bash
FILE=outc1i.txt
C=0
THRDS="1 4 8 12 16"
RFNMT="3 4 5 6 7"
NRUNS=1

usage(){
    echo "Usage: $0 [number of test runs]"
}

header(){
    # Write ARGS directive
    echo ARGS Refinement Threads
    # Write caption depending on block size
	if [ $C = "1" ]
	then
	    echo CAPTION:2-by-2 block structure | tee -a ${FILE}
	else
	    echo CAPTION:3-by-3 block structure | tee -a ${FILE}
	fi
}

run(){

    for refine in $RFNMT
    do
        for ompt in $THRDS
        do
             echo "Arg Refinements $refine Threads $ompt"
             export OMP_NUM_THREADS=$ompt
             ./elastic -r $refine -c $C
        done
    done
}

if [[ $# -eq 1 ]]; then
        NRUNS=$1
fi

header | tee -a ${FILE}

for rns in $(seq 1 $NRUNS)
do
        echo "+++++++++========== RUN $rns ==========+++++++++"
        run | tee -a ${FILE}
done