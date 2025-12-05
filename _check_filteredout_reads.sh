#!/bin/bash


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage="Usage: $0 -d <output_dir> -s <bam_file>"

while getopts ":d:s:" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    s)
      bam_file="$OPTARG"
      ;;
    \?)
      echo -e "Invalid option: -$OPTARG. \n Usage: $usage" >&2
      exit 1
      ;;
  esac
done

filtered_bam_file=${output_dir}/entropy_filtered.bam

thrownaway_reads_bam="${output_dir}/entropy_thrownaway_reads.bam"

# (Negative control)
# Check reads that correspond to thrown-away barcodes (filtered out by entropy threshold)
# 
echo "Check filtered-out reads in BAM file ....... "
time { samtools view $bam_file | cut -f 1 | sort | uniq > ${output_dir}/all_readnames.txt ; }
time { samtools view $filtered_bam_file | cut -f 1 | sort | uniq > ${output_dir}/entropy_filtered_readnames.txt ; }

# # Get the read names that were filtered out
time { comm -23 ${output_dir}/all_readnames.txt ${output_dir}/entropy_filtered_readnames.txt > ${output_dir}/filteredout_readnames.txt ; }

# Extract the reads with the filtered-out read names into a separate BAM file
samtools view -N ${output_dir}/filteredout_readnames.txt $bam_file -bo ${thrownaway_reads_bam}


echo "Check filtered-out reads in BAM file ....... Done"