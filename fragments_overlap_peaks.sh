#!/bin/bash


while getopts ":d:" opt; do
  case $opt in
    d) output_dir="$OPTARG" 
      ;;
    \?) 
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;  esac
done

comparison_subdir=${output_dir}/peaks_comparison
lost_peaks_file=$comparison_subdir/lost_peaks.bed

peaks_subdir=${output_dir}/peaks
peaks_file=${peaks_subdir}/peaks.bed

entropy_peaks_subdir=${output_dir}/peaks_entropy_filtered
entropy_peaks_file=${entropy_peaks_subdir}/peaks.bed



frag_file=${output_dir}/fragments.tsv


frag_peak_overlap_subdir=${output_dir}/fragments_overlap_peaks
mkdir -p $frag_peak_overlap_subdir

frag_overlap_lostpeaks_file=$frag_peak_overlap_subdir/overlap_lost_peaks.bed


bedtools intersect -a $frag_file -b $lost_peaks_file  -u > $frag_overlap_lostpeaks_file
# count the number of fragments overlapping with lost peaks per each 
cut -f4 $frag_overlap_lostpeaks_file | sort | uniq -c > $frag_peak_overlap_subdir/overlap_lost_peaks_counts.txt 

bedtools intersect -a $frag_file -b $peaks_file  -u > $frag_peak_overlap_subdir/overlap_peaks.bed
cut -f4 $frag_peak_overlap_subdir/overlap_peaks.bed | sort | uniq -c > $frag_peak_overlap_subdir/overlap_peaks_counts.txt

bedtools intersect -a $frag_file -b $entropy_peaks_file  -u > $frag_peak_overlap_subdir/overlap_entropy_peaks.bed
cut -f4 $frag_peak_overlap_subdir/overlap_entropy_peaks.bed | sort | uniq -c > $frag_peak_overlap_subdir/overlap_entropy_peaks_counts.txt


# inplace trailing off the leading spaces in the output files
sed -i 's/^[ ]*//'  $frag_peak_overlap_subdir/overlap_lost_peaks_counts.txt  
sed -i 's/^[ ]*//'  $frag_peak_overlap_subdir/overlap_peaks_counts.txt
sed -i 's/^[ ]*//'  $frag_peak_overlap_subdir/overlap_entropy_peaks_counts.txt