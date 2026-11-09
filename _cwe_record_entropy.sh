#!/bin/bash

set -euo pipefail

_DefaultNprocess=4 # Default number of processes for entropy calculation, if not specified by user
_WindowSize=3000 # Default window size for fragmentation, if not specified by user
THIS_DIR=$(dirname $(realpath ${BASH_SOURCE[0]}))

# echo ${THIS_DIR}
# exit 0
source "${THIS_DIR}/config.sh"
source ${PYTHON_ENV}

Usage="Usage: $0 -o output_dir -g species [-w window_size (default: ${_WindowSize})]"

output_dir=""
window_size=""
nprocess=""
species=""

while getopts "o:w:g:p:" opt; do
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
    p)
      nprocess="$OPTARG"
      ;;
    *)
      echo "${Usage}"  && exit 1
      ;;
  esac
done

if [[ -z ${window_size} ]]; then
    window_size=${_WindowSize}
fi

if [[ -z ${nprocess} ]]; then
    nprocess=${_DefaultNprocess}
fi

if [[ -z ${species} ]]; then
    echo "Species must be specified with -g option." >&2
    echo "${Usage}" >&2 &&    exit 1
else
    species_genome_size_file="${THIS_DIR}/ref/${species}_genome_chromsize.tsv"
    blacklist_file="${THIS_DIR}/ref/${species}/${species}-blacklist.bed"
    if [[ ! -f ${species_genome_size_file} ]]; then
        echo "Genome size file for species ${species} not found: ${species_genome_size_file}" >&2 && exit 1
    fi
fi

if [[ -z ${output_dir} ]]; then
    echo "Output directory must be specified with -o option." >&2
    echo "${Usage}" >&2 && exit 1
fi

if [[ ! -d ${output_dir}/chromosome_fragments ]]; then
    echo "No chromosome_fragments/ in ${output_dir}. Run sade.sh on this output directory first." >&2 && exit 1
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
    if [[  -f ${output_dir}/chromosome_fragments/${chromosome}_fragments.tsv ]]; then
        if [[ ! -f "${output_dir}/chromosome_fragments/.${chromosome}_insertion_done.log" ]]; then
            python "${THIS_DIR}/count_perchrom_tn5_insertions.py" \
                --output_dir="${output_dir}" \
                --genome_chromsize_file="${species_genome_size_file}" \
                --chromosome="${chromosome}" \
                --blacklist_file="${blacklist_file}" \
                --windowsize="${window_size}" \
                --shortcut_mode "yes"  && touch "${output_dir}/chromosome_fragments/.${chromosome}_insertion_done.log"
        else
            echo "Insertion for chromosome ${chromosome} already counted. Skip to calculate entropy."
        fi

      # Calculate entropy for the specified chromosome
      python "${THIS_DIR}/calculate_entropy.py" --output_dir="${output_dir}" --chromosome="${chromosome}"  --genomesize_file="${species_genome_size_file}" --shortcut_mode "yes" --nprocess="${nprocess}"
    fi
done