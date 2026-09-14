source ~/python_virtuenv/atac_tobias3.13/bin/activate


output_dir="/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison"

frag_bam_file=$output_dir/sorted_tagged_union_bc_fragments.bam

res_dir=${output_dir}/cmp_SC/entropy_peak_atac_cluster
rna_cluster_file=$res_dir/cmp_SC_entropypeakatac_cluster.csv

clusters="0 1 2 3 4 5 6 7 8 9 10 11"
clusters="11"
clusters="0 1 2 3 4 5 6 7 8 9 10"
for c in $clusters;
do 
	awk -F '\t' -v OFS='\t' -v clus="$c" '{ if ($5==clus) print $1 }' $rna_cluster_file > $res_dir/atac_c_${c}.txt
	samtools view -h $frag_bam_file | awk -F '\t' -v OFS='\t'  '
     NR==FNR {
	bc_list[$1]; next
      }

	{ if ($1 ~ /^@/) { print $0; next }
	if (  $1 in bc_list ) print $0  
	}
	'  $res_dir/atac_c_${c}.txt  - |  samtools sort -@ 8 -o $res_dir/sorted_farg_atac_c_${c}.bam 
	
	samtools index $res_dir/sorted_farg_atac_c_${c}.bam
	bamCoverage -b $res_dir/sorted_farg_atac_c_${c}.bam \
		    -o $res_dir/sorted_farg_atac_c_${c}.bw \
		    -bs 1 \
		    -p 8

done
