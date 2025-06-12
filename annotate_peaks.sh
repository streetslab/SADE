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
lost_peaks=$comparison_subdir/lost_peaks.bed


annotation_subdir="${output_dir}/annotation"
mkdir -p $annotation_subdir
annotated_newly_discovered_peaks=$annotation_subdir/annotated_newly_discovered_peaks.bed
annotated_lost_peaks=$annotation_subdir/annotated_lost_peaks.bed


# annotate peaks differences 
annotatePeaks.pl $newly_discovered_peaks $genome_name -annStats  $annotation_subdir/"annStats_newly_discovered_peaks.txt" > $annotated_newly_discovered_peaks
annotatePeaks.pl $lost_peaks $genome_name -annStats  $annotation_subdir/"annStats_lost_peaks.txt" > $annotated_lost_peaks



# Annotate peaks 
peak_calling_subdir="${output_dir}/peaks"
peaks_bed_file="${peak_calling_subdir}/peaks.bed"

entropy_peak_calling_subdir="${output_dir}/peaks_entropy_filtered"
entropy_peaks_bed_file="${entropy_peak_calling_subdir}/peaks.bed"

annotatePeaks.pl $peaks_bed_file $genome_name -annStats $annotation_subdir/"annStats_peaks.txt" > $annotation_subdir/annotated_peaks.bed
annotatePeaks.pl $entropy_peaks_bed_file $genome_name -annStats $annotation_subdir/"annStats_entropy_peaks.txt" > $annotation_subdir/annotated_entropy_peaks.bed