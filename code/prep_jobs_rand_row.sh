#!/bin/bash

source ./package_params_rand_row.sh

cd $FOLDER
python prep_jobs_rand_row.py -d $DATAFILE -d2 $RANDDATAFILE -a $ARGSFILE -t $TEST -o $OUTPUTNAME -n $JOBNUM -p $PARALLELNUM -f $FDR -c $COEFNUM

echo
echo Done prepping jobs