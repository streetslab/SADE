#!/bin/bash

_WindowSize=3000 # Default window size for fragmentation, if not specified by user
THIS_DIR=$(dirname $(realpath ${BASH_SOURCE[0]}))

# echo ${THIS_DIR}
# exit 0
source "${THIS_DIR}/config.sh"
source ${PYTHON_ENV}

while getopts "o:w:g:" opt; do
  case ${opt} in
    o)
      output_dir=$OPTARG
      output_dir=${output_dir%/} # Remove trailing slash if present
      ;;
    w)
      window_size=$OPTARG
      ;;
    g)
      species="$OPTARG"
      ;;
    *)
      echo "Usage: $0 -o output_dir  -w window_size -g species"  && exit 1
      ;;
  esac
done

if [[ -z ${window_size} ]]; then
    window_size=${_WindowSize}
fi

if [[ -z ${species} ]]; then
    echo "Species must be specified with -g option."
    exit 1
else
    species_genome_size_file="${THIS_DIR}/ref/${species}_genome_chromsize.tsv"
    if [[ ! -f ${species_genome_size_file} ]]; then
        echo "Genome size file for species ${species} not found: ${species_genome_size_file}"
        exit 1
    fi
fi


chromosome_list=()
while IFS=$'\t' read -r chrom _res; do
    if [[ $chrom =~ ^'#' ]]; then
        continue
    else
        chromosome_list+=("$chrom")
    fi
done < "${species_genome_size_file}"




for chromosome in "${chromosome_list[@]}"; do
    echo "Processing chromosome: ${chromosome}"
    if [[ -f ${output_dir}/chromosome_fragments/${chromosome}_fragments.tsv ]]; then
      python "${THIS_DIR}/count_perchrom_tn5_insertions.py" \
          --output_dir="${output_dir}" \
          --genome_chromsize_file="${species_genome_size_file}" \
          --chromosome="${chromosome}" \
          --windowsize="${window_size}" \
          --shortcut_mode "yes"

      # Calculate entropy for the specified chromosome
      python "${THIS_DIR}/calculate_entropy.py" --output_dir="${output_dir}" --chromosome="${chromosome}"  --genomesize_file="${species_genome_size_file}" --shortcut_mode "yes"
    fi
done