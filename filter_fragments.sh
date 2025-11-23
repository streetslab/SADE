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
entropy_df_file="${output_dir}/${chromosome}_barcode_entropy_df.tsv"

# File to save results
filtered_bc_df_file="${output_dir}/entropy_filtered_bc_df.tsv"
filtered_fragments_file="${output_dir}/filtered_fragments.tsv"



entropy_threshold=$(awk -F',' 'NR==1 {print $2}' "${entropy_cutoff_file}")


# Only keep barcodes with entropy >= threshold
echo "Filtering fragments corresponding to barcodes passed the entropy threshold............"

temp_bc_file="${output_dir}/bc_pass_entropy.tsv"
awk -F',' -v threshold="${entropy_threshold}" 'NR==1 || $2 >= threshold'  ${entropy_df_file} > ${filtered_bc_df_file}
awk -F',' 'NR>1 {print $1}' ${filtered_bc_df_file} > ${temp_bc_file}
grep -f ${temp_bc_file} ${fragment_file} > ${filtered_fragments_file}
awk -v OFS='' -v prefix='CB:Z:' '{{print prefix, $1}}' ${temp_bc_file} > ${output_dir}/temp_bc_CBZ.txt

echo "Done filtering fragments."


