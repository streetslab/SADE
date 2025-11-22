#!/bin/bash

while getopts "f:o:" opt; do
  case $opt in
    f) 
      fragment_file="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    *)
      echo "Usage: $0 -f <fragment_file> -o <output_dir>" >&2
      exit 1
      ;;
  esac
done

# make output directory if it doesn't exist
mkdir -p $output_dir

awk '{print $4}' $fragment_file | sort | uniq > ${output_dir}/unique_barcodes.tsv