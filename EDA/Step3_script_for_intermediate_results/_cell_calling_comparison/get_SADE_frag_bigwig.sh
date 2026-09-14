source ~/python_virtuenv/atac_tobias3.13/bin/activate


output_dir="/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison"

frag_bam_file=$output_dir/sorted_tagged_union_bc_fragments.bam

res_dir=${output_dir}/Entropy_cells_peaks
SADE_bc_file=$res_dir/Entropy_filtered_bc_CBZ.txt

	samtools view -h $frag_bam_file | awk -F '\t' -v OFS='\t'  '
     NR==FNR {
	bc_list[$1]; next
      }

	{ if ($1 ~ /^@/) { print $0; next }
	if (  $1i2 in bc_list ) print $0  
	}
	'  ${SADE_bc_file}   - |  samtools sort -@ 8 -o $res_dir/SADE_bc_fragments.bam 
	
	samtools index $res_dir/SADE_bc_fragments.bam  
	bamCoverage -b $res_dir/SADE_bc_fragments.bam \
		    -o $res_dir/SADE_bc_fragments.bw \
		    -bs 1 \
		    -p 18

