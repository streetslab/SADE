#!/bin/bash



while getopts "d:" opt; do
  case $opt in
    d) 
      output_dir="$OPTARG"
        output_dir=${peak_dir%/}  # remove trailing slash if exists
      ;;
      ;;
    \?) 
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


comparison_subdir="${output_dir}/peaks_comparison"
mkdir -p $comparison_subdir


# those sub_dir were pre-difined in process.sh
entropy_peak_calling_subdir="${output_dir}/peaks_entropy_filtered"
peak_calling_subdir="${output_dir}/peaks"


peaks_file=$peak_calling_subdir/peaks.bed
peaks_entropy_file=$entropy_peak_calling_subdir/peaks.bed



bedtools intersect -wo -a $peaks_file -b $peaks_entropy_file  >  $comparison_subdir/peaks_intersect.bed