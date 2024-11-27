#!/usr/bin/env bash

piscem_executable=$1 ## rust executable for piscem rust executable
map_rad=$2 ## rad file output by piscem
whitelist_file=$3 ## whitelist file
rev_comp=$4
threads=$5
output_path=$6

echo "Correcting barcodes"
/usr/bin/time -o $output_path/time_whitelist.txt $piscem_executable generate-permit-list \
    --input $map_rad \
    --unfiltered-pl $whitelist_file \
    --rev-comp $rev_comp \
    --output-dir $output_path

echo "Sort"
/usr/bin/time -o $output_path/time_sort.txt $piscem_executable sort \
    --input-dir $output_path \
    --rad-dir $map_rad \
    --threads $threads
