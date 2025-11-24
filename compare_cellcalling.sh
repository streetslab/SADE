#!/bin/bash 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPT_DIR}/config.sh



while getopts ":d:b:c:" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    b)
      crbarcode_file="$OPTARG"
      ;;
    c)
      chromosome="$OPTARG" # has default value
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


# Plot entropy values of barcodes overlaying with CellRanger cell-calling labels
source ${PYTHON_ENV}
python ${SCRIPT_DIR}/overlay_entropy_CRcelllabel_plot.py \
    --output_dir $output_dir \
    --chromosome $chromosome \
    --crbarcode_file $crbarcode_file 