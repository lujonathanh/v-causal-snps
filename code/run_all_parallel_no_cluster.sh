while read script; do
    echo Submitting Parallel Script $script
    sh -c "$script"
done <parallel_script_list.txt