#!/bin/bash


while getopts "s:o:" opt; do
  case $opt in
    s)
      bam_file="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


mkdir -p $output_dir

output_peak_file="${output_dir}/peaks_w_blacklistregion.bed"
output_bedgraph_file="${output_dir}/bedgraph.bed"
output_bed_file="${output_dir}/intervals.bed"

echo "Call peaks with Genrich .... "  && \
echo " ... could take a while, ..." && \
sorted_bam="${output_dir}/sorted_bam.bam"  && \
samtools sort -n  -@ 6 $bam_file -o $sorted_bam  # use 6 threads for sorting
#TODO: Optimize on number of threads used for samtools sort


$Genrich -t $sorted_bam -o $output_peak_file  \
        -k $output_bedgraph_file \
        -b $output_bed_file \
        -r \
        -j  && rm $sorted_bam \
        &&  echo "Call peaks with Genrich .... Done"