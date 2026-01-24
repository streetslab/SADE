#!/bin/bash

output_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000_mixD'
chromosome='chr1'


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


entropy_cutoff_file="${output_dir}/entropy_cutoff.csv"
fragment_file="${output_dir}/fragments.tsv"
chromosome_fragment_dir=${output_dir}/chromosome_fragments
entropy_df_file="${output_dir}/${chromosome}_barcode_entropy_df.tsv"

# File to save results
filtered_bc_df_file="${output_dir}/entropy_filtered_bc_df.tsv"
filtered_fragments_file="${output_dir}/filtered_fragments.tsv"



entropy_threshold=$(awk -F',' 'NR==1 {print $2}' "${entropy_cutoff_file}")


# Only keep barcodes with entropy >= threshold
echo "Filtering fragments corresponding to barcodes passed the entropy threshold............"

temp_bc_file="${output_dir}/bc_pass_entropy.tsv"
awk -F',' -v threshold="${entropy_threshold}" 'NR==1 || $2 >= threshold'  ${entropy_df_file} > ${filtered_bc_df_file}
awk -F',' 'NR>1 {print $1}' ${filtered_bc_df_file} > ${temp_bc_file}  # Remove header for grep step 
# ###### Approach 1 grep 
# grep -f ${temp_bc_file} ${fragment_file} > ${filtered_fragments_file}  # Single threaded version #TODO: update command options to have both versions
# ###### Approach 2: parallel grep
# #part 1: Create directory to hold per-chromosome fragment files
# filtered_perchrom_frag_dir="${output_dir}/filtered_chromosome_fragments"
# if [ ! -d ${filtered_perchrom_frag_dir} ] ; then
#     mkdir -p ${filtered_perchrom_frag_dir}
# fi
# #part 2: parallel grep
# Parallelize on chromosome_fragments files for speedup. -- Grep is still slow 
# find ${chromosome_fragment_dir} -name "*_fragments.tsv" -print0 | xargs -I{}  -0 -P 4 \
#     sh -c 'res_file=$(basename $3) ; LC_ALL=C  grep -f $1 $3 > "${2}/${res_file}" ' -- \
#     "${temp_bc_file}" "${filtered_perchrom_frag_dir}" {}
# cat ${filtered_perchrom_frag_dir}/*_fragments.tsv > ${filtered_fragments_file}
# Clean up temporary files
# rm -rf ${filtered_perchrom_frag_dir}

# ###### Approach 3: use rg (ripgrep) -- fastest
rg -f ${temp_bc_file}  ${fragment_file} > ${filtered_fragments_file}



awk -v OFS='' -v prefix='CB:Z:' '{print prefix, $1}' "${temp_bc_file}" > "${output_dir}/Entropy_filtered_bc_CBZ.txt"

echo "Done filtering fragments."


