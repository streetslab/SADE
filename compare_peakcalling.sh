#!/bin/bash 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh


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
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done




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
if [ ! -f ${entropy_peak_calling_subdir}/peaks.bed ] ; then
    bash ${SCRIPT_DIR}/call_peaks.sh \
    -s ${filtered_bam_file} \
    -o ${entropy_peak_calling_subdir} 
fi

# MODULE 4:  (II) peak calling on the original BAM file
peak_calling_subdir="${output_dir}/peaks"
if [ ! -f ${peak_calling_subdir}/peaks.bed ] ; then
    bash ${SCRIPT_DIR}/call_peaks.sh \
    -s ${bam_file} \
    -o ${peak_calling_subdir}
fi


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

