#!/bin/bash

source ./package_params.sh

cd $FOLDER
python prep_jobs.py -d $DATAFILE -a $ARGSFILE -t $TEST -o $OUTPUTNAME -n $JOBNUM

echo
echo Done prepping jobs