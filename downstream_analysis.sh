#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh


while getopts "o:b:" opt; do
  case $opt in
    o) output_dir="$OPTARG" 
       output_dir=${output_dir%/} # Remove trailing slash if present
       ;; 
    b) barcode_file="$OPTARG" 
        ;;
    *) echo "Usage: $0 -o <output_dir> -b <barcode_file>" >&2; exit 1 ;;
  esac
done



filtered_fragments_file="${output_dir}/filtered_fragments.tsv"

DownstreamReanalysis_dir="${output_dir}/DownstreamReanalysis"
mkdir -p ${DownstreamReanalysis_dir}



## Step 1: Sort filtered fragments in .gz format and create a sorted and indexed file in .tbi 
filtered_fragments_sorted_gz_file="${DownstreamReanalysis_dir}/filtered_fragments_sorted.tsv.gz"
if [[ ! -f ${filtered_fragments_sorted_gz_file} ]]; then
    # Sort the filtered fragments file and keep only the first three columns (chromosome, start, end)
    sort -k1,1 -k2,2n  ${filtered_fragments_file} > ${DownstreamReanalysis_dir}/filtered_fragments_sorted.tsv
    # Compress to .gz format -- this is in-place compression 
    bgzip ${DownstreamReanalysis_dir}/filtered_fragments_sorted.tsv  
    # Create an indexed file for the sorted fragments
    tabix -p bed ${filtered_fragments_sorted_gz_file}
fi


## Step 2: Run get_pic.r to generate the entropy peak count matrix
# Activate the conda environment that has R 
DownstreamReanalysis_dir="${output_dir}/DownstreamReanalysis"
entropy_bc_peak_count_h5="${DownstreamReanalysis_dir}/entropy_bc_peak_count.h5"
if [[ ! -f ${entropy_bc_peak_count_h5} ]]; then
    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate ${Renv_Conda} 
    # Run the R script to get the PIC matrix 
    Rscript ${SCRIPT_DIR}/get_pic.r --barcode_file="${barcode_file}" --res_dir="${output_dir}"

    conda deactivate
fi

## Step 3: Analyze the PIC matrix to compare with RNA's clustering/cell-type annotation's consistency

