#!/bin/bash 

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

comparison_subdir="${output_dir}/peaks_comparison"
newly_discovered_peaks=$comparison_subdir/newly_discovered_peaks.bed 


annotation_subdir="${output_dir}/annotation"
mkdir -p $annotation_subdir
annotated_newly_discovered_peaks=$annotation_subdir/annotated_newly_discovered_peaks.bed


annotatePeaks.pl $newly_discovered_peaks $genome_name > $annotated_newly_discovered_peaks