#!/bin/bash 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


while getopts ":d:g:" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    g)
      genome_name="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


# MODULE 1. Compare peaks called between before and post entropy filtering
comparison_subdir=${output_dir}/peaks_comparison
newly_discovered_peaks=$comparison_subdir/newly_discovered_peaks.bed
if [ ! -f $newly_discovered_peaks ] ; then
    bash ${SCRIPT_DIR}/compare_peaks.sh -d $output_dir
fi



# MODULE 2. Annotate peaks with functional regions 
bash ${SCRIPT_DIR}/annotate_peaks.sh -d $output_dir -g $genome_name