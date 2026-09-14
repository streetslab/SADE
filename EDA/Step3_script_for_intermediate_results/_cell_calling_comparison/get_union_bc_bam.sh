
awk -F '\t' -v OFS="" -v prefix=CB:Z: 'NR>1 {print prefix, $1, -1}' Union_cell_3set.tsv  > Union_cells_bc_CBZ.txt
#awk -F '\t' -v OFS='\t' '{print $0, $1"-1"}' Union_cell_3set_with_RNA_cluster_and_DNAdebrisflag.tsv > _Union_cell_3set_with_RNA_cluster_and_DNAdebrisflag.tsv 
# ^ keep a column for barocde with '-1' appendix for filtering fragment file

bamfile='/mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/VIB_10xmultiome_2/outs/possorted_bam.bam'
pres_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'

union_dir=${pres_dir}/_cell_calling_comparison

bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
                        -o ${union_dir} \
                        -f ${union_dir}/Union_cells_bc_CBZ.txt
