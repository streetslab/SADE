#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ${SCRIPT_DIR}/config.sh


Msg="Usage: $0 -o output_dir -c chromosome -w window_size -g species"

while getopts "o:c:w:g:" opt; do
  case ${opt} in
    o)
      output_dir=$OPTARG
      output_dir=${output_dir%/} # Remove trailing slash if present
      ;;
    c)
      chromosome=$OPTARG
      ;;
    w)
      window_size=$OPTARG
      ;;
    g)
      species="$OPTARG"
      ;;
    *)
      echo $Msg  && exit 1
      ;;
  esac
done


# Step 1: Split fragments file by chromosome 
chromosome_fragment_dir=${output_dir}/chromosome_fragments
finish_file=${chromosome_fragment_dir}/.split_fragments_by_chromosome.finished
if [[ ! -f ${finish_file} ]]; then
    echo "Splitting fragments file by chromosome..."
    bash $SCRIPT_DIR/split_fragments_by_chrom.sh -o ${output_dir}  && \
    touch ${finish_file} # create finish file
fi


# Step 2: Count Tn5 insertion frequency for each cell barcode
## This can be done in parallel for multiple chromosomes in python
## TODO: add parallelization 
source ${PYTHON_ENV}
insert_frequency_file="${output_dir}/{chromosome}_insert_frequency.pickle"
species_genome_size_file="${species}_genome_chromsize.tsv"
if [[ ! -f ${insert_frequency_file} ]]; then
    python ${SCRIPT_DIR}/count_perchrom_tn5_insertions.py \
        --output_dir="${output_dir}" \
        --genome_chromosize_file="${SCRIPT_DIR}/ref/${species_genome_size_file}" \
        --chromosome="${chromosome}" \
        --windowsize="${window_size}"
fi


# Step 3: Calculate entropy for each chromosome separately
source ${PYTHON_ENV}
barcode_entropy_df_file="${output_dir}/${chromosome}_barcode_entropy_df.tsv"
if [[ ! -f ${barcode_entropy_df_file} ]]; then
    # Calculate entropy for the specified chromosome
    python ${SCRIPT_DIR}/calculate_entropy.py --output_dir="${output_dir}" --chromosome="${chromosome}"
fi

