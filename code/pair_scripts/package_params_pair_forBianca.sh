#!/bin/bash
export OUTPUTNAME=edgeR-reg-avg-3-randomized

export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-avg/featurecounts.genes.TPM.selected_reps.ln.surrogate_variables_corrected.protein_coding-edgeR-reg-randomized-avg.txt
export DATAFILE2=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-avg/featurecounts.genes.TPM.selected_reps.ln.surrogate_variables_corrected.protein_coding-edgeR-reg.txt

export ARGSFILE=granger_args-3-noproc.txt

export TEST=gp
export JOBNUM=60
export PARALLELNUM=15
export TIMEUPPERBOUND=8
export MEMUPPERBOUND=15000
export FOLDER=della/$OUTPUTNAME

# Above is in minutes


echo FOLDER is $FOLDER
echo DATAFILE is $DATAFILE
echo DATAFILE2 is $DATAFILE2
echo ARGSFILE is $ARGSFILE
echo TEST is $TEST
echo OUTPUTNAME is $OUTPUTNAME
echo JOBNUM is $JOBNUM
echo PARALLELNUM is $PARALLELNUM
echo TIMEUPPERBOUND is $TIMEUPPERBOUND
echo MEMUPPERBOUND is $MEMUPPERBOUND