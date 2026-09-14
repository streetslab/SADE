

output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000F'
res_dir="${output_dir}/_cell_calling_comparison"
frag_file=$res_dir/union_bc_fragments.bed


cmp_list2="$res_dir/Rawpeaks $res_dir/Entropy_cells_peaks $res_dir/FRIP_cells_peaks $res_dir/TSS_cells_peaks"
# cmp_list2="$res_dir/Union_cells_peaks"
# cmp_list2="$res_dir/Commoncells_peaks"
for comparison_subdir in $cmp_list2
do 
	peaks=$comparison_subdir/peaks.bed
	frag_out_dir=${comparison_subdir}/frag_overlap
	mkdir -p ${frag_out_dir}
	bedtools intersect -a $frag_file -b $peaks -u > $frag_out_dir/frag_in_peaks.tsv

done



