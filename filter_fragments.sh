#!/bin/bash

Msg="Usage: $0 -o output_dir  -s genome_saturation_cutoff"

while getopts "o:s:" opt; do
  case ${opt} in
    o)
      output_dir=$OPTARG
      output_dir=${output_dir%/} # Remove trailing slash if present
      ;;
    s)
      genome_saturation_cutoff=$OPTARG
      ;;
    *)
      echo "$Msg" && exit 1
      ;;
  esac
done

if [[ -z ${output_dir} ||  -z ${genome_saturation_cutoff} ]]; then
    echo "$Msg" && exit 1
fi

fragment_file="${output_dir}/fragments.tsv"
entropy_cutoff_file="${output_dir}/entropy_cutoff.csv"
entropy_df_file="${output_dir}/calculated_barcode_entropy_df.tsv"

# File to save results
filtered_bc_wflag_df_file="${output_dir}/entropy_filtered_bc_w_DNAdebrisflag_df.tsv"
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


entropy_rank_cutoff=$(awk -F',' 'NR==3 {print $2}' "${entropy_cutoff_file}")
echo $entropy_rank_cutoff

# Only keep barcodes with entropy >= threshold
echo "Filtering fragments corresponding to barcodes passed the entropy threshold............"

temp_bc_file="${output_dir}/bc_pass_entropy.tsv"
tmp_df_file="${output_dir}/tmp.tsv"
# S1. Filter barcodes passed entropy threshold. 
awk -F',' -v OFS=',' 'NR==1 {print $0}'  "${entropy_df_file}" > "${tmp_df_file}" # Write header
awk -F ','  -v OFS=','  'NR>1 {print $0}'  "${entropy_df_file}"   | sort -k2,2 -g -r -t ',' | head -n "${entropy_rank_cutoff}"  >> "${tmp_df_file}"
# S2a. Flag barcodes "YES" for DNA debris by estimated-P_closed_state < (1 - genome_saturation_cutoff)
awk -F ',' -v OFS=',' -v cutoff="${genome_saturation_cutoff}" \
    'NR==1 {print $0, "DNA_debris"}  NR>1 {if  ($7 < (1-cutoff))  print $0, "YES" ; else print $0, "NO" }'  "${tmp_df_file}" > "${filtered_bc_wflag_df_file}"
# S2b. Final filtered barcodes: passed entropy threshold and not flagged as DNA debris 
awk -F',' -v OFS=',' 'NR==1 {print $0}'  "${filtered_bc_wflag_df_file}" > "${filtered_bc_df_file}"   # Write header
awk -F ',' -v OFS=',' '$8=="NO" {print $0}'  "${filtered_bc_wflag_df_file}" >> "${filtered_bc_df_file}"
# S3. Add DNA debris flag to the original entropy_df_file for all barcodes (including those not passed entropy threshold)
awk -F ',' -v OFS=',' 'NR==FNR {dna_debris[$1]=$8; next} NR==1 {print $0, "DNA_debris"} NR>1 {flag=(dna_debris[$1]!="" ? dna_debris[$1] : "NO"); print $0, flag}'  "${filtered_bc_wflag_df_file}" "${entropy_df_file}" > "${output_dir}/calculated_barcode_entropy_wDNAdebrisflag_df.tsv"


awk -F',' 'NR>1 {print $1}' "${filtered_bc_df_file}" > "${temp_bc_file}"  # Remove header for grep step 


# ###### Approach 3: use rg (ripgrep) -- fastest
rg -f "${temp_bc_file}"  "${fragment_file}" > "${filtered_fragments_file}"


awk -v OFS='' -v prefix='CB:Z:' '{print prefix, $1}' "${temp_bc_file}" > "${output_dir}/Entropy_filtered_bc_CBZ.txt"

# Clean up temporary files
rm "${tmp_df_file}"
echo "Done filtering fragments."


