while read script; do
    echo Submitting Script $script
    ./$script
done <script_list.txt