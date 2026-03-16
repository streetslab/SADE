#!/bin/bash


while getopts "d:" opt; do
  case $opt in
    d) output_dir="$OPTARG" 
       output_dir=${output_dir%/} # Remove trailing slash if present
       ;;
    *) echo "Usage: $0 -o <output_dir> -b <barcode_file>" >&2; exit 1 ;;
  esac
done


DownstreamReanalysis_dir="${output_dir}/DownstreamReanalysis"

filtered_fragments_file="${output_dir}/filtered_fragments.tsv"

# Sort the filtered fragments file and keep only the first three columns (chromosome, start, end)
sort -k1,1 -k2,2n  ${filtered_fragments_file} > ${DownstreamReanalysis_dir}/filtered_fragments_sorted.tsv
# Compress to .gz format -- this is in-place compression 
bgzip ${DownstreamReanalysis_dir}/filtered_fragments_sorted.tsv  
# Create an indexed file for the sorted fragments in .tbi 
tabix -p bed ${DownstreamReanalysis_dir}/filtered_fragments_sorted.tsv.gz