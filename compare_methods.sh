#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

while getopts ":d:b:s:g:f:n" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    b)
      crbarcode_file="$OPTARG"
      ;;
    s)
      bam_file="$OPTARG"
      ;;
    g)
      genome_name="$OPTARG"
      ;;
    f)
      fragments_file_gz="$OPTARG"
      ;;
    n)
      sample_name="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done



## CR filtering results
mkdir -p "${output_dir}/_CR_FRIP"
cp $crbarcode_file "${output_dir}/_CR_FRIP/barcodes.tsv"  # copy as barcodes.tsv 


## TSS filtering results 
mkdir -p "${output_dir}/_ArchR_TSS"
conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
source "${conda_dir}/etc/profile.d/conda.sh" 
conda activate ${Renv_Conda} 
# Run the R script to get ArchR filtering results 
Rscript ${SCRIPT_DIR}/get_tss_w_archr.r -f="${fragments_file_gz}" --res_dir="${output_dir}" --sample_name="${sample_name}" --genome_name="${genome_name}"

conda deactivate



## identify raw peaks using all the reads in the experiment
calling_cmp_dir="${output_dir}/_cell_calling_comparison"
Raw_peak_dir="${calling_cmp_dir}/Rawpeaks" 
mkdir -p "${Raw_peak_dir}"
bash ${SCRIPT_DIR}/call_peaks.sh \
  -s ${bam_file} \
  -o ${Raw_peak_dir}
bedtools subtract -a "${Raw_peak_dir}/peaks_w_blacklistregion.bed"  -b "${SCRIPT_DIR}/ref/${genome_name}/${genome_name}-blacklist.bed"  > "${Raw_peak_dir}/peaks.bed"
