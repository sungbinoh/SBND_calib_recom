#!/bin/bash

python_script="make_root_list_with_dir.py"
base_dir="/data/sungbino/sbnd/2024_data/calib/nov_dec"
output_prefix="list_dune_gpu01"

runs=(
    "17625"
    "17632"
    "17635"
    "17647"
    "17648"
    "17657"
    "17658"
    "17659"
    "17660"
    "17661"
    "17662"
    "17663"
    "17664"
    "17665"
)

for run in "${runs[@]}"; do
    input_dir="${base_dir}/run_${run}"
    output_file="${output_prefix}_run_${run}.txt"
    
    echo "Processing: $input_dir -> $output_file"
    python "$python_script" "$input_dir" "$output_file"
done

echo "All runs processed."
