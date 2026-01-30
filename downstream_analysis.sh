#!/bin/bash
set -e # exit on error and dont continue
set -o pipefail # catch errors in piped commands

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh

Usage="Usage: $0 -d <output_dir> -b <crbarcode_file>
          -d <output_dir: Directory where the output files are located from running auto_process.sh. Required>
          -b <crbarcode_file: CellRanger cell barcode file (barcodes.tsv) to use for comparison. Required>
          "


while getopts "d:b:" opt; do
  case $opt in
    d) output_dir="$OPTARG" 
       output_dir=${output_dir%/} # Remove trailing slash if present
       ;; 
    b) barcode_file="$OPTARG" 
        ;;
    *) echo "${Usage}" >&2 && exit 1 ;;
  esac
done


if [[ -z ${output_dir} ]]; then
    echo -e "Please provide output directory with -d \n"
    echo "${Usage}" >&2  && exit 1
fi

if [[ -z ${barcode_file} ]]; then
    echo -e "Please provide CellRanger barcode file with -b \n"
    echo "${Usage}" >&2  && exit 1
fi



DownstreamReanalysis_dir="${output_dir}/DownstreamReanalysis"
mkdir -p ${DownstreamReanalysis_dir}



## Step 1: Sort filtered fragments in .gz format and create a sorted and indexed file in .tbi 
filtered_fragments_sorted_gz_file="${DownstreamReanalysis_dir}/filtered_fragments_sorted.tsv.gz"
if [[ ! -f ${filtered_fragments_sorted_gz_file} ]]; then
    bash ${SCRIPT_DIR}/sort_filterd_frag.sh -d ${output_dir} 
fi


## Step 2: Run get_pic.r to generate the entropy peak count matrix
# Activate the conda environment that has R 
entropy_bc_peak_count_h5="${DownstreamReanalysis_dir}/entropy_bc_peak_count.h5ad"
if [[ ! -f ${entropy_bc_peak_count_h5} ]]; then
    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate ${Renv_Conda} 
    # Run the R script to get the PIC matrix 
    Rscript ${SCRIPT_DIR}/get_pic.r --barcode_file="${barcode_file}" --res_dir="${output_dir}"

    conda deactivate
fi

## Step 3: Analyze the PIC matrix using PEAKVI
source ${PYTHON_ENV_PEAKVI}
checkpoint_file="${DownstreamReanalysis_dir}/downstream_complete.log"
if [[ ! -f $checkpoint_file ]]; then
    # Run the PEAKVI analysis script
    python ${SCRIPT_DIR}/peackvi_atac.py --res_dir="${output_dir}"
fi

