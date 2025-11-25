#!/bin/bash 

chromosome='chr1'
# << GLOBAL VARIABLES >>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh


while getopts ":d:" opt; do
  case $opt in
    d) output_dir="$OPTARG" 
       output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    \?) 
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;  
  esac
done


# MODULE 1. Compare fragments overlapping% with called peaks 
frag_peak_overlap_subdir=${output_dir}/fragments_overlap_peaks
frag_overlap_entropy_peaks_file=$frag_peak_overlap_subdir/overlap_entropy_peaks_counts.txt
if [ ! -f $frag_overlap_entropy_peaks_file ] ; then
    echo "Running fragments overlap with peaks analysis..."
    bash ${SCRIPT_DIR}/fragments_overlap_peaks.sh -d $output_dir 
fi 


# MODULE 2. Make plots of fragments overlapping with peaks
source ${PYTHON_ENV}
python ${SCRIPT_DIR}/visualize_fragments_overlap_peaks.py --res_dir $output_dir --chromosome $chromosome






#### TODO: ARCHIVE this file.. ALL The modules are in compare_peakcalling.sh now.  REMOVE this file later.  ####