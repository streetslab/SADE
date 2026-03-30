set -x

output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000F'


PYTHON_ENV='/home/syyang/python_virtuenv/atac_tobias3.13/bin/activate'
source $PYTHON_ENV

rfilter='/home/syyang/GitRepo/atac/ref/hg38/hg38-blacklist.bed'

amulet_path='/home/syyang/GitRepo/AMULET'
amulet="${amulet_path}/AMULET.sh"


des_dir="${output_dir}/_cell_calling_comparison"
Union_BC_file="${des_dir}/Union_cells_bc.txt"
filtered_fragments_file="${des_dir}/union_bc_fragments.tsv.gz"

Union_BC_csv="${des_dir}/union_filtered_bc.csv"
touch "${Union_BC_file}"  # Create empty file if it doesn't exist
echo 'barcode,is__cell_barcode' > "${Union_BC_csv}"
awk -v OFS=',' '{print $1, 1}' "${Union_BC_file}" >> "${Union_BC_csv}"


amulet_output=${des_dir}/AMULET_output
mkdir -p "${amulet_output}"




$amulet  "${filtered_fragments_file}"  "${Union_BC_csv}"  \
    ${amulet_path}/human_autosomes.txt \
"${rfilter}"  "${amulet_output}"  "${amulet_path}"
