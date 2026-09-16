

output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000F'
res_dir="${output_dir}/_cell_calling_comparison"
frag_file=$res_dir/union_bc_fragments.bed

blacklist_file='/home/syyang/GitRepo/SADE/ref/hg38/hg38-blacklist.bed'


bedtools intersect -a $frag_file -b $blacklist_file -u > $res_dir/frag_in_blacklist.tsv


