

output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000F'
output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000_test_norm'

filtered_bc_file=${output_dir}/Entropy_filtered_bc_CBZ.txt
fragment_file="${output_dir}/fragments.tsv"

data_dir='/mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/VIB_10xmultiome_2'
CR_BC_file="${data_dir}/outs/filtered_peak_bc_matrix/barcodes.tsv"



des_dir="${output_dir}/_UnionBC_fragments"
mkdir -p "${des_dir}"   

Entropy_BC_file="${des_dir}/Entropy_filtered_bc.txt"
Union_BC_file="${des_dir}/UnionBC.txt"
Union_BC_csv="${des_dir}/UnionBC.csv"
filtered_fragments_file="${des_dir}/filtered_fragments.tsv"

awk  -F ':'  '{print $3}' "${filtered_bc_file}" > "${Entropy_BC_file}"


# # Record union Barcodes from CR and Entropy cell calling results
cat "${Entropy_BC_file}" "${CR_BC_file}" | sort | uniq > "${Union_BC_file}"

# filter fragments corresponding to the union barcodes
rg -f "${Union_BC_file}"  "${fragment_file}" > "${filtered_fragments_file}"
gzip "${filtered_fragments_file}"


