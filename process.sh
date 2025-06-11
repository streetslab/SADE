#!/bin/bash 

chromosome='chr1'
# << GLOBAL VARIABLES >>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# environment setup
source ${SCRIPT_DIR}/config.sh
# load $PYTHON_ENV
#echo "python_path $PYTHON_ENV" # for debugging


## BEFORE ANYTHING ELSE: process options 
while  getopts "f:o:b:" opt; do
  case $opt in
    f) 
      fragment_file="$OPTARG"
      ;;
    o)
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


# make sure to only execute the script on the fragment file that's within the output directory
frag_file=${output_dir}/$(basename $fragment_file)
if [[  $fragment_file == *.gz ]]; then
    frag_file=${frag_file%.gz} # update frag_file name
fi

if [ ! -f $frag_file ] ; then
    mkdir -p $output_dir
    gunzip -c $fragment_file > ${frag_file}
fi


# MODULE 1:  Calculate entropies of each barcode
source ${PYTHON_ENV}

entropy_file="${output_dir}/${chromosome}_barcode_entropy.pickle"
if [ ! -f ${entropy_file} ]; then
    python ${SCRIPT_DIR}/calculate_entropy.py \
        --res_dir $output_dir \
        --frag_file $frag_file \
        --genome_chromsize "${SCRIPT_DIR}/ref/human_genome_chromsize.tsv" --chromosome "${chromosome}"

fi

# MODULE 2:  Plot entropy of each barcode and overlay with CR cell-calling labels
python ${SCRIPT_DIR}/overlay_entropy_CRcelllabel_plot.py \
        --res_dir $output_dir \
        --frag_file $frag_file \
        --entropy_file $entropy_file \
        --crbarcode_file $crbarcode_file
