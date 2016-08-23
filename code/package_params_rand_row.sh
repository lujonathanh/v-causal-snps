#!/bin/bash



export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/small_TPM-avg.txt
#/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-avg/featurecounts.hg.TPM.selected_reps.ln.surrogate_variables_corrected.protein_coding-edgeR-reg-avg.txt
export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/small_TPM-randomized-avg.txt
#/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-avg/featurecounts.hg.TPM.selected_reps.ln.surrogate_variables_corrected.protein_coding-edgeR-reg-randomized-avg.txt

export GENES=small_TPM-fdr-test
#prot
export DEG=
#er
export SAMPLE=
#avg


# These should correspond! PGC = g, Enet = e
export ARGSFILE=enet_args-2.txt
export CAUSAL=enet
export TEST=e
export COEFNUM=2
export FDR="0.05"


export OUTPUTNAME="$GENES-$DEG-$SAMPLE-$CAUSAL-$COEFNUM"

export NNODES=5
export PPN=16
export PARALLELNUM=$(($NNODES * $PPN))


# only used for deciding how many jobs to prep, not in actual program
export NROWS=100
#3545
export NPERSCRIPT=7
export JOBNUM=$(expr $NROWS / $NPERSCRIPT)

export TIMEUPPERBOUND=61
export MEMUPPERBOUND=45000
export FOLDER=della/$OUTPUTNAME
# Above is in minutes











echo FOLDER is $FOLDER
echo DATAFILE is $DATAFILE
echo RANDDATAFILE is $RANDDATAFILE
echo ARGSFILE is $ARGSFILE
echo TEST is $TEST
echo OUTPUTNAME is $OUTPUTNAME
echo JOBNUM is $JOBNUM
echo PARALLELNUM is $PARALLELNUM
echo TIMEUPPERBOUND is $TIMEUPPERBOUND
echo MEMUPPERBOUND is $MEMUPPERBOUND