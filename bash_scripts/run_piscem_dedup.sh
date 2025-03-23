#!/usr/bin/env bash

af_exec=$1 ## rust executable for piscem rust executable
map_rad=$2 ## directory containing rad file with alignments output by piscem
whitelist_file=$3 ## whitelist file
rev_comp=$4
threads=$5
output_path=$6
i=$7

echo "Correcting barcodes"
/usr/bin/time -o $output_path/time_whitelist_threads=${threads}_${i}.out $af_exec generate-permit-list \
    --input $map_rad \
    --unfiltered-pl $whitelist_file \
    --rev-comp $rev_comp \
    --output-dir $output_path

echo "Sorting"
/usr/bin/time -o $output_path/time_sort_threads=${threads}_${i}.out $af_exec sort \
    --input-dir $output_path \
    --rad-dir $map_rad \
    --threads $threads
