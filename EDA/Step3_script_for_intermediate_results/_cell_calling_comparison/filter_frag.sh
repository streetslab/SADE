
frag_file='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/fragments.tsv'
# For union bcs 

awk -F '\t' -v OFS=""  'NR>1 {print $1}' Union_cell_3set_all_info.tsv  > Union_cells_bc.txt


rg -f  Union_cells_bc.txt $frag_file > union_bc_fragments.tsv

sort -k1,1 -k2,2n union_bc_fragments.tsv  > union_bc_fragments.bed

#awk -F '\t' -v OFS='\t' '{print $0, $1"-1"}' Union_cell_3set_with_RNA_cluster_and_DNAdebrisflag.tsv > _Union_cell_3set_with_RNA_cluster_and_DNAdebrisflag.tsv 
# ^ keep a column for barocde with '-1' appendix for filtering fragment file
