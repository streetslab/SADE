#!/bin/bash 
set -e # exit on error and dont continue
set -o pipefail # catch errors in piped commands

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh

Usage="Usage: $0 -d <output_dir> -g <genome_name> -s <bam_file>
  -d: <output_dir: Directory where the output files are located from running auto_process.sh. Required>
  -g: <genome_used_for_read_mapping_that_resulted_fragments_file. Required: 'hg38', 'mm10' etc.> 
  -s: <Original BAM file used for peak calling. Required>
  "

while getopts ":d:g:s:" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    g)
      genome_name="$OPTARG"
      ;;
    s)
      bam_file="$OPTARG"
      ;;
    \?)
      echo -e "Invalid option: -$OPTARG. \n Usage: $Usage" >&2
      exit 1
      ;;
  esac
done

if [[ -z ${output_dir} || -z ${genome_name} || -z ${bam_file} ]]; then
    echo "Missing required arguments."
    echo -e "$Usage" >&2
    exit 1
fi

if [[ -z ${bam_file} ]]; then
    echo "Please provide original BAM file used for peak calling with -s"
    echo -e "$Usage" >&2
    exit 1
fi


parameter_file="${output_dir}/parameters.csv"

# extract chromosome from parameter file
chromosome_line=$(grep "chromosome" ${parameter_file})
chromosome=${chromosome_line#*,}


# MODULE 3:  Filter BAM file with barcodes passed the entropy threshold
filtered_bc_file=${output_dir}/Entropy_filtered_bc_CBZ.txt  #check-point for MODULE 2
filtered_bam_file=${output_dir}/entropy_filtered.bam #check-point for MODULE 3
if [ ! -f ${filtered_bam_file} ] ; then
    bash ${SCRIPT_DIR}/filter_bam_with_barcodes.sh \
        -s ${bam_file} \
        -o ${output_dir} \
        -f ${filtered_bc_file}
fi


# MODULE 4:  (I) peak calling on the filtered BAM file
entropy_peak_calling_subdir="${output_dir}/peaks_entropy_filtered"
if [ ! -f ${entropy_peak_calling_subdir}/peaks_w_blacklistregion.bed ] ; then
    bash ${SCRIPT_DIR}/call_peaks.sh \
    -s ${filtered_bam_file} \
    -o ${entropy_peak_calling_subdir} 
fi

# MODULE 4:  (II) peak calling on the original BAM file
peak_calling_subdir="${output_dir}/peaks"
if [ ! -f ${peak_calling_subdir}/peaks_w_blacklistregion.bed ] ; then
    bash ${SCRIPT_DIR}/call_peaks.sh \
    -s ${bam_file} \
    -o ${peak_calling_subdir}
fi


# MODULE 4: Remove blacklist regions from called peaks
bash ${SCRIPT_DIR}/remove_blacklist_region.sh \
    -o ${output_dir} \
    -g ${genome_name}



# Now ---  Compare the peak calling results


# MODULE 5. Compare peaks called between before and post entropy filtering
comparison_subdir=${output_dir}/peaks_comparison
newly_discovered_peaks=$comparison_subdir/newly_discovered_peaks.bed
if [ ! -f $newly_discovered_peaks ] ; then
    bash ${SCRIPT_DIR}/compare_peaks.sh -d $output_dir
fi



# MODULE 6. Annotate peaks with functional regions 
annotation_subdir=${output_dir}/annotation
annStats_peaks=$annotation_subdir/annStats_peaks.txt
if [ ! -f $annStats_peaks ] ; then
    bash ${SCRIPT_DIR}/annotate_peaks.sh -d $output_dir -g $genome_name
fi


# MODULE 7. Make plots of peaks comparisons
source ${PYTHON_ENV}
python ${SCRIPT_DIR}/visualize_annotation.py --output_dir $output_dir



# MODULE 8. Compare fragments overlapping% with called peaks 
frag_peak_overlap_subdir=${output_dir}/fragments_overlap_peaks
frag_overlap_entropy_peaks_file=$frag_peak_overlap_subdir/overlap_entropy_peaks_counts.txt
if [ ! -f $frag_overlap_entropy_peaks_file ] ; then
    echo "Running fragments overlap with peaks analysis..."
    bash ${SCRIPT_DIR}/fragments_overlap_peaks.sh -d $output_dir 
fi 


# MODULE 9. Make plots of fragments overlapping with peaks
source ${PYTHON_ENV}
python ${SCRIPT_DIR}/visualize_fragments_overlap_peaks.py --res_dir $output_dir --chromosome $chromosome

