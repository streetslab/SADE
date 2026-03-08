#!/bin/bash

Msg="Usage: $0 -o output_dir -c chromosome"

while getopts "o:c:" opt; do
  case ${opt} in
    o)
      output_dir=$OPTARG
      output_dir=${output_dir%/} # Remove trailing slash if present
      ;;
    c)
      chromosome=$OPTARG
      ;;
    *)
      echo "$Msg" && exit 1
      ;;
  esac
done

if [[ -z ${output_dir} || -z ${chromosome} ]]; then
    echo "$Msg" && exit 1
fi

fragment_file="${output_dir}/fragments.tsv"
entropy_cutoff_file="${output_dir}/entropy_cutoff.csv"
entropy_df_file="${output_dir}/calculated_barcode_entropy_df.tsv"

# File to save results
filtered_bc_df_file="${output_dir}/entropy_filtered_bc_df.tsv"
filtered_fragments_file="${output_dir}/filtered_fragments.tsv"



# Validate required input files exist
for file in "${fragment_file}" "${entropy_cutoff_file}" "${entropy_df_file}"; do
    if [[ ! -f "${file}" ]]; then
        echo "ERROR: Required file not found: ${file}" >&2
        echo "Please ensure upstream pipeline steps completed successfully." >&2
        exit 1
    fi
done

# entropy_threshold=$(awk -F',' 'NR==1 {print $2}' "${entropy_cutoff_file}")
entropy_rank_cutoff=$(awk -F',' 'NR==3 {print $2}' "${entropy_cutoff_file}")
echo $entropy_rank_cutoff

# Only keep barcodes with entropy >= threshold
echo "Filtering fragments corresponding to barcodes passed the entropy threshold............"

temp_bc_file="${output_dir}/bc_pass_entropy.tsv"
awk -F',' -v OFS=',' 'NR==1 {print $1, $2}'  "${entropy_df_file}" > "${filtered_bc_df_file}"  # Write header
awk -F ','  -v OFS=','  'NR>1 {print $1, $2}'  "${entropy_df_file}"   | sort -k2,2 -r -t ',' | head -n "${entropy_rank_cutoff}"  >> "${filtered_bc_df_file}"  
awk -F',' 'NR>1 {print $1}' "${filtered_bc_df_file}" > "${temp_bc_file}"  # Remove header for grep step 


# ###### Approach 3: use rg (ripgrep) -- fastest
rg -f "${temp_bc_file}"  "${fragment_file}" > "${filtered_fragments_file}"



awk -v OFS='' -v prefix='CB:Z:' '{print prefix, $1}' "${temp_bc_file}" > "${output_dir}/Entropy_filtered_bc_CBZ.txt"

echo "Done filtering fragments."


