source ~/python_virtuenv/atac_tobias3.13/bin/activate


output_dir="/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison"

frag_bam_file=$output_dir/sorted_tagged_union_bc_fragments.bam

pres_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'


res_dir=${output_dir}/TSS_bc_lowFRIP_bc/individual_bc_frag

res_dir=${output_dir}/common_bc_highFRIP_bc/individual_bc_frag
bc="AAACTCGCAATCCTAC"
bc_index=1

	bc_list_file=$res_dir/$f
	bw_dir=${bc_list_file/%.tsv/}
	mkdir -p ${bw_dir}

		bc1=${bc}"-1"


		samtools view -h $frag_bam_file | awk -F '\t' -v OFS='\t' -v BC="$bc1" '
		{ if ($1 ~ /^@/) { print $0; next }
		if (  $1 == BC ) print $0  
		}
		'  - |  samtools sort -@ 8 -o $bw_dir/sorted_farg_bc_${bc_index}.bam 
	
		samtools index $bw_dir/sorted_farg_bc_${bc_index}.bam
		bamCoverage -b $bw_dir/sorted_farg_bc_${bc_index}.bam \
		    -o $bw_dir/sorted_farg_bc_${bc_index}.bw \
		    -bs 1 \
		    -p 16
