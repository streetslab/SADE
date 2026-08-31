#!/bin/bash 
set -o errexit # exit on error and dont continue
set -o pipefail # catch errors in piped commands

# GLOBAL VARIABLES >>>
_CHROMOSOME='chr1'
_WindowSize=3000 # default window size
_CandidateInflectionPoints=2    # default number of candidate inflection points to consider for entropy threshold selection.
_GenomeSaturationCutoff=0.5 #  _GenomeSaturationCutoff := (1 - estimated_closed_region_genome_coverage). Higher value, less strict filtering of DNA-debris. 
                            # Human cells typically have <0.7 genome coverage in scATAC-seq. this cutoff is pretty loose. 
# GLOBAL VARIABLES <<<

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# environment setup
source "${SCRIPT_DIR}/config.sh"


Usage="Usage: $0 
          -f <scATACseq fragments_file. Required. This file has '\t' as field separator (.tsv or .bed) and has at least 4 columns "chr", "start-pos", "end-pos", and "droplet-barcode".> 
          -o <output_dir. Required> 
          -g <genome_used_for_read_mapping_that_resulted_fragments_file. Required: 'hg38', 'mm10' etc.> 
          [-c <chromosome>. Default: ${_CHROMOSOME}] 
          [-w <window_size>. Default: ${_WindowSize}]
          [-k <candidate_inflection_points>. Default: ${_CandidateInflectionPoints}]
          [-s <genome_saturation_cutoff>. Default: ${_GenomeSaturationCutoff}. Values above cutoff means the genome is too saturated with fragments, 
                                          they correspond to DNA debris.
                                          Mathematically valid value is in the range of (0, 1], but do not specify manually unless you know the 
                                          average genome coverage in basepair (bp) of a random given cell across all cell types for the species in the sample.
                                          Due to this reason, we manually set range to be [0.1, 1], where 1 means no filtering of DNA debris.]
          "

## BEFORE ANYTHING ELSE: process options 
while  getopts "f:o:g:c:w:k:s:" opt; do
  case $opt in
    f) 
      fragment_file="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    g)
      species="$OPTARG"
      # currently: 'hg38', 'mm10'
      ;;
    c)
      chromosome="$OPTARG" # has default value
      ;;
    w)
      window_size="$OPTARG"  # has default value
      ;;
    k)
      candidate_inflection_points="$OPTARG" # has default value
      ;;
    s)
      genome_saturation_cutoff="$OPTARG" # has default value
      ;;
    \?) 
      echo -e "Invalid option: -$OPTARG \n$usage" >&2
      exit 1
      ;;
  esac
done


# set default value if not provided
## set chromosome
if [[ -z ${chromosome} ]]; then
    chromosome=${_CHROMOSOME}
fi
## size window size 
if [[ -z ${window_size} ]]; then
    window_size=${_WindowSize}
fi
## candidate inflection points for entropy threshold selection
if [[ -z ${candidate_inflection_points} ]]; then
    candidate_inflection_points=${_CandidateInflectionPoints}
fi
## 
if [[ -z ${fragment_file} ]]; then
    echo "Please provide fragments file with -f"
    echo "${Usage}" >&2  && exit 1
fi
##
if [[ -z ${species} ]]; then
    echo "Please provide species with -g"
    echo "${Usage}" >&2  && exit 1
fi
##
if [[ -z ${output_dir} ]]; then
    echo "Please provide directory to save results"
    echo "${Usage}" >&2  && exit 1
fi
##
if [[ -z ${genome_saturation_cutoff} ]]; then
    genome_saturation_cutoff=${_GenomeSaturationCutoff}
fi

# Environment validation >>>
# Validate Python environment
if [[ -z "${PYTHON_ENV}" ]]; then
    echo "ERROR: PYTHON_ENV not set in config.sh"
    echo "Please edit config.sh and set PYTHON_ENV to your Python virtual environment activate script"
    exit 1
fi

if [[ ! -f "${PYTHON_ENV}" ]]; then
    echo "ERROR: Python environment file not found: ${PYTHON_ENV}"
    echo "Please check PYTHON_ENV setting in config.sh"
    exit 1
fi
# Environment validation <<<

# Parameter validation >>>
cutoff=''
cutoff=$(awk -v specified_cutoff="${genome_saturation_cutoff}"  -v default_cutoff="${_GenomeSaturationCutoff}" 'BEGIN { 
    if (specified_cutoff < 0.1 || specified_cutoff > 1) {
        print default_cutoff 
    }
}')
if [[ ! -z ${cutoff} ]]; then
    genome_saturation_cutoff=${_GenomeSaturationCutoff}
fi
# Parameter validation <<<


figures_subdir="${output_dir}/figures"
mkdir -p "${figures_subdir}"

# make sure to only execute the script on the fragment file that's within the output directory
frag_file="${output_dir}/fragments.tsv"
checkpoint_file="${output_dir}/_finish_format_fragments.log"
if [[ ! -f ${checkpoint_file} ]] ; then
    echo "Checking fragments file format for entropy calculation..."

    if [[ ${fragment_file} == *.gz ]]; then
        cat "${fragment_file}" | zcat | awk -F '\t' '
            BEGIN { valid=0; invalid=0 }
            {
                if (NF == 5) {
                    print $0
                    valid++
                } else {
                    invalid++
                }
            }
            END {
                if (valid == 0) {
                    print "ERROR: No valid 5-column rows found in fragments file" > "/dev/stderr"
                    exit 1
                }
                if (invalid > 0) {
                    print "WARNING: Filtered out " invalid " invalid rows (" valid " valid rows kept)" > "/dev/stderr"
                }
            }' - > ${frag_file} && touch ${checkpoint_file} 
    elif [[ ${fragment_file} == *.tsv ]]; then 
        # sometimes fragments.tsv from CR does not have the correct record.. 
        awk -F '\t' '
            BEGIN { valid=0; invalid=0 }
            {
                if (NF == 5) {
                    print $0
                    valid++
                } else {
                    invalid++
                }
            }
            END {
                if (valid == 0) {
                    print "ERROR: No valid 5-column rows found in fragments file" > "/dev/stderr"
                    exit 1
                }
                if (invalid > 0) {
                    print "WARNING: Filtered out " invalid " invalid rows (" valid " valid rows kept)" > "/dev/stderr"
                }
            }' "${fragment_file}" > "${frag_file}"  &&  touch "${checkpoint_file}"  # create empty file as checkpoint
    else
        echo "Unsupported fragment file format. Please provide .tsv or .tsv.gz file." >&2 && exit 1
    fi
fi



# Record parameters used in this run
echo "Parameters used in this run:" > "${output_dir}/parameters.csv"
echo "species,${species}" >> "${output_dir}/parameters.csv"
echo "chromosome,${chromosome}" >> "${output_dir}/parameters.csv"
echo "window_size,${window_size}" >> "${output_dir}/parameters.csv"
echo "genome_saturation_cutoff,${genome_saturation_cutoff}" >> "${output_dir}/parameters.csv"
echo "candidate_inflection_points,${candidate_inflection_points}" >> "${output_dir}/parameters.csv"

# MODULE 1:  Calculate entropies of each barcode
barcode_entropy_df_file="${output_dir}/calculated_barcode_entropy_df.tsv"
if [ ! -f ${barcode_entropy_df_file} ] ; then
  bash ${SCRIPT_DIR}/calculate_entropy.sh \
          -o "${output_dir}" \
          -g "${species}" \
          -c "${chromosome}" \
          -w "${window_size}"
fi 


# ++MODULE 1.5: automatically determine entropy threshold based on knee plot
auto_entropy_file="${output_dir}/entropy_cutoff.csv"
if [ ! -f ${auto_entropy_file} ] ; then
    source ${PYTHON_ENV}
    python ${SCRIPT_DIR}/find_threshold.py \
        --output_dir "${output_dir}" \
        --chromosome "${chromosome}" \
        --k  "${candidate_inflection_points}"

    # Record entropy threshold into the parameters file
    awk 'NR==1' ${auto_entropy_file} >> ${output_dir}/parameters.csv
fi



# MODULE 2:  Filter fragments based on entropy thresholdded barcodes
#            and plot it overlaying with CR cell-calling labels
filtered_frag_file="${output_dir}/filtered_fragments.tsv"  #check-point for MODULE 2
filtered_bc_file="${output_dir}/Entropy_filtered_bc_CBZ.txt"  #check-point for MODULE 2
if [ ! -f ${filtered_frag_file} ] ; then 
    # Filter fragments based on the entropy threshold
    bash "${SCRIPT_DIR}/filter_fragments.sh" \
        -o "${output_dir}" \
        -s "${genome_saturation_cutoff}"
fi 


# Visualization 
source ${PYTHON_ENV}
python "${SCRIPT_DIR}/refine_figures.py" \
    --output_dir "${output_dir}" \
    --chromosome "${chromosome}"

