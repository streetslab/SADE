#!/bin/bash 

## BEFORE ANYTHING ELSE: process options 
while  getopts ":d:" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    \?) 
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

annotation_subdir=${output_dir}/annotation
entropy_peaks=$annotation_subdir/annotated_entropy_peaks.bed
newly_discovered_peaks=$annotation_subdir/annotated_newly_discovered_peaks.bed
lost_peaks=$annotation_subdir/annotated_lost_peaks.bed

igv_dir=${annotation_subdir}/igv
mkdir -p $igv_dir

awk -F '\t' -v OFS='\t' '{print $2, $3, $4, $1, $16, $8}' $newly_discovered_peaks > $igv_dir/annotated_newly_discovered_peaks.bed
awk -F '\t' -v OFS='\t' '{print $2, $3, $4, $1, $16, $8}' $lost_peaks > $igv_dir/annotated_lost_peaks.bed
awk -F '\t' -v OFS='\t' '{print $2, $3, $4, $1, $16, $8}' $entropy_peaks > $igv_dir/annotated_entropy_peaks.bed

grep 'promoter-TSS' $igv_dir/annotated_newly_discovered_peaks.bed  > $igv_dir/_temp.bed
awk -F 'S ' '{print $1"S"}' $igv_dir/_temp.bed  > $igv_dir/annotated_newly_discovered_TSS_peaks.bed
rm $igv_dir/_temp.bed

