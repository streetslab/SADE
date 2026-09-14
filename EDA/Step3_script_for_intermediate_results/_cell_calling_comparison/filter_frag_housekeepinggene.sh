
res_dir="/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison"

frag_file="$res_dir/union_bc_fragments.bed"
housekeeping_gene_file="$res_dir/_housekeepinggenes/pbmc_housekeeping_genes.bed"

frag_out_dir=$res_dir

bedtools intersect -a $frag_file -b $housekeeping_gene_file -wo > $frag_out_dir/frag_in_housekeeping.bed
