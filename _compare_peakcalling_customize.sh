#!/bin/bash 
set -e # exit on error and dont continue
set -o pipefail # catch errors in piped commands

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh

Usage="Usage: $0 -d <output_dir> -g <genome_name> [-s <bam_file>]
  -d: <output_dir: Directory where the output files are located from running auto_process.sh. Required>
  -g: <genome_used_for_read_mapping_that_resulted_fragments_file. Required: 'hg38', 'mm10' etc.>
  -s: <Original BAM file. Required for PEAK_CALLER=genrich, not needed for PEAK_CALLER=macs> 
  -b: <barcode_file: filtered cell barcode file (barcodes.tsv) to use for comparison. Required>
  -r: <result_dir>. Relative path that will be appended to output_dir to store the comparison results. Default: 'comparison_results'>
  "

while getopts ":d:g:s:b:r:" opt; do
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
    b)
      cmpbarcode_file="$OPTARG"
      ;;
    r)
      result_subdir="$OPTARG"
      result_subdir=${result_subdir%/}  # remove trailing slash if exists
      result_dir="${output_dir}/${result_subdir}"
      ;;
    \?)
      echo -e "Invalid option: -$OPTARG. \n Usage: $Usage" >&2
      exit 1
      ;;
  esac
done

if [[ -z ${output_dir} || -z ${genome_name} ]]; then
    echo "Missing required arguments."
    echo -e "$Usage" >&2
    exit 1
fi

# BAM file is only required for genrich peak caller
if [[ "${PEAK_CALLER}" != "macs" && -z ${bam_file} ]]; then
    echo "Please provide original BAM file used for peak calling with -s (required for PEAK_CALLER=${PEAK_CALLER:-genrich})"
    echo -e "$Usage" >&2
    exit 1
fi



mkdir -p ${result_dir}

entropy_peak_calling_subdir="${result_dir}/peaks_entropy_filtered"
othermethod_peak_calling_subdir="${result_dir}/peaks"


if [[ "${PEAK_CALLER}" == "macs" ]]; then
    # Not implemented -- pointless. MACS has way too many false positive peaks. 
    echo 'skip'
else
    # Genrich mode: requires BAM file + samtools (original behavior)

    # MODULE 3:  Filter BAM file with specified barcodes file
    filtered_bc_file=${result_dir}/othermethod_filtered_bc_CBZ.txt  #check-point for MODULE 2

    awk -v OFS='' -v prefix='CB:Z:' '{print prefix, $1}' "${cmpbarcode_file}" > "${filtered_bc_file}"

    filtered_bam_file=${result_dir}/othermethod_filtered.bam #check-point for MODULE 3
    filtered_sam_file=${result_dir}/othermethod_filtered.sam # intermediate file for MODULE 3
    if [ ! -f ${filtered_bam_file} ] ; then

        # Save the header lines
        samtools view -H $bam_file > $filtered_sam_file
        # Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
        echo "  Filtering alignments in BAM file ............"
        time { samtools view $bam_file | LC_ALL=C rg -j 0 -F -f $filtered_bc_file >> $filtered_sam_file ;  } 
        time { samtools view -b $filtered_sam_file > $filtered_bam_file ; }
        rm $filtered_sam_file
    fi

    # MODULE 4:  (I) peak calling on the SADE entropy filtered BAM file
    # Already done in standard pipeline. Just need to copy the results here for comparison.
    cp -r "${output_dir}/peaks_entropy_filtered" "${result_dir}/peaks_entropy_filtered"

    # MODULE 4:  (II) peak calling on the comparison cellc-calling method filtered BAM file
    if [ ! -f ${othermethod_peak_calling_subdir}/peaks_w_blacklistregion.bed ] ; then
        bash ${SCRIPT_DIR}/call_peaks.sh \
        -s ${filtered_bam_file} \
        -o ${othermethod_peak_calling_subdir}
    fi
fi


# MODULE 4: Remove blacklist regions from called peaks
bash ${SCRIPT_DIR}/remove_blacklist_region.sh \
    -o ${result_dir} \
    -g ${genome_name}



# Now ---  Compare the peak calling results


# MODULE 5. Compare peaks called between before and post entropy filtering
comparison_subdir=${result_dir}/peaks_comparison
newly_discovered_peaks=$comparison_subdir/newly_discovered_peaks.bed
if [ ! -f $newly_discovered_peaks ] ; then
    bash ${SCRIPT_DIR}/compare_peaks.sh -d $result_dir
fi



# MODULE 6. Annotate peaks with functional regions 
annotation_subdir=${result_dir}/annotation
annStats_peaks=$annotation_subdir/annStats_peaks.txt
if [ ! -f $annStats_peaks ] ; then
    bash ${SCRIPT_DIR}/annotate_peaks.sh -d $result_dir -g $genome_name
fi


# MODULE 7. Make plots of peaks comparisons
source ${PYTHON_ENV}
python ${SCRIPT_DIR}/visualize_annotation.py --output_dir $result_dir



# # MODULE 8. Compare fragments overlapping% with called peaks 
# frag_peak_overlap_subdir=${result_dir}/fragments_overlap_peaks
# frag_overlap_entropy_peaks_file=$frag_peak_overlap_subdir/overlap_entropy_peaks_counts.txt
# if [ ! -f $frag_overlap_entropy_peaks_file ] ; then
#     echo "Running fragments overlap with peaks analysis..."
#     bash ${SCRIPT_DIR}/fragments_overlap_peaks.sh -d $result_dir 
# fi 


