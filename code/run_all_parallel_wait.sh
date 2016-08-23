while read script; do
    echo Submitting Parallel Script $script
    clusterize -l $TIMEUPPERBOUND:00 -m $MEMUPPERBOUND -n $NNODES -p $PPN  -c "time $script & wait"
done <parallel_script_list.txt