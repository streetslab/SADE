#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


while getopts "o:g:" opt; do
  case $opt in
    o)
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


peak_calling_subdir="${output_dir}/peaks"
entropy_peak_calling_subdir="${output_dir}/peaks_entropy_filtered"

# Remove blacklist regions from called peaks. 
# Doing this step post peak calling so that we do not have to adjust p-values in Genrich. 
# See: https://github.com/jsh58/Genrich/issues/23
bedtools subtract -a ${entropy_peak_calling_subdir}/peaks_w_blacklistregion.bed  -b ${SCRIPT_DIR}/ref/${genome_name}/${genome_name}-blacklist.bed  > ${entropy_peak_calling_subdir}/peaks.bed
bedtools subtract -a ${peak_calling_subdir}/peaks_w_blacklistregion.bed  -b ${SCRIPT_DIR}/ref/${genome_name}/${genome_name}-blacklist.bed  > ${peak_calling_subdir}/peaks.bed