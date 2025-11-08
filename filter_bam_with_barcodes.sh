#!/bin/bash


while getopts "s:o:f:" opt; do
  case $opt in
    s)
      BAM_FILE="$OPTARG"
      ;;
    o)
      OUTPUT_DIR="$OPTARG"
      ;;
    f)
      FILTER_BC_FILE="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


echo "Filter bam file with barcodes that passed the entropy threshold ----> "
echo " ... could take a while, ..."


filtered_sam_file="${OUTPUT_DIR}/entropy_filtered.sam"
filtered_bam_file="${OUTPUT_DIR}/entropy_filtered.bam"

# Save the header lines
samtools view -H $BAM_FILE > $filtered_sam_file

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
# samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILTER_BC_FILE >> $filtered_sam_file
time { samtools view $BAM_FILE | LC_ALL=C rg -j 0 -F -f $FILTER_BC_FILE >> $filtered_sam_file ;  } # multi-threaded grep using ripgrep
# Convert filtered.sam to BAM format
samtools view -b $filtered_sam_file > $filtered_bam_file


rm $filtered_sam_file


echo "Filter bam file with barcodes that passed the entropy threshold ---->  Done"