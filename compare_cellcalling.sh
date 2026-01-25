#!/bin/bash 
set -e # exit on error and dont continue
set -o pipefail # catch errors in piped commands

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh


Usage="Usage: $0 -d output_dir -b crbarcode_file [-c chromosome]
  -d <output_dir: Directory where the output files are located from running auto_process.sh. Required>
  -b <crbarcode_file: CellRanger cell barcode file (barcodes.tsv) to use for comparison. Required>
"

while getopts ":d:b:" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    b)
      crbarcode_file="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [[ -z ${output_dir} ]]; then
    echo -e "Please provide output directory with -d \n"
    echo "${Usage}" >&2  && exit 1
fi

if [[ -z ${crbarcode_file} ]]; then
    echo -e "Please provide CellRanger barcode file with -b \n"
    echo "${Usage}" >&2  && exit 1
fi

parameter_file="${output_dir}/parameters.csv"

# extract chromosome from parameter file
chromosome_line=$(grep "chromosome" ${parameter_file})
chromosome=${chromosome_line#*,}

# MODULE 0:  PREPEARATION for global variables and files to plot
## Remove # lines from the fragments file and count the number of unique fragments per barocode
frag_file=${output_dir}/fragments.tsv
checkpoint_file="${output_dir}/.finish_module_0"
if [[ ! -f ${checkpoint_file} ]] ; then
  sed '/^#/d' ${frag_file}| cut -f4 | sort | uniq -c > ${output_dir}/total_fragments_counts.txt
  ## Inplace trailing off the leading spaces in the output files
  sed -i 's/^[ ]*//'  ${output_dir}/total_fragments_counts.txt 
  touch ${checkpoint_file}  # create empty file as checkpoint for MODULE 0
fi



# Plot entropy values of barcodes overlaying with CellRanger cell-calling labels
source ${PYTHON_ENV}
echo "Make figures comparing Entropy cell-calling with CellRanger cell-calling .... "  && \
python ${SCRIPT_DIR}/overlay_entropy_CRcelllabel_plot.py \
    --output_dir $output_dir \
    --chromosome $chromosome \
    --crbarcode_file $crbarcode_file 