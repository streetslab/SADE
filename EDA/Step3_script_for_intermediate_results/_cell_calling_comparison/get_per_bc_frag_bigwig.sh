source ~/python_virtuenv/atac_tobias3.13/bin/activate


output_dir="/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison"

frag_bam_file=$output_dir/sorted_tagged_union_bc_fragments.bam


# For specific barcode
frag_bam_dir="${output_dir}/per_bc_frag_bam2"
common_bc_list_select="AAACAAGCAAGTTAGT AAACAAGCATTTCTTC AAACCGCCAACGGTAC AAACGAACACCTGTAT AAACTAGTCGGCTATC AAACTCGCAATCCTAC AAAGCCTAGGATTATG AAAGCCTAGTTACGTT AAAGGCAAGGCTTGTC AAAGGCGCACTAATCC"
for bc in $common_bc_list_select
do
        f="${bc}-1_fragments.bam"
        if [[ ! -f "${frag_bam_dir}/${f/%bam/bw}" ]]
        then
        bamCoverage -b  "${frag_bam_dir}/$f"  \
                    -o "${frag_bam_dir}/${f/%bam/bw}"  \
                    -bs 1 \
                    -p 18
        fi
done

exit 0


# For specific set -- obsolete
res_dir=${output_dir}/cmp_SC/individual_bc_frag

clusters="DNA_debris_bc_list.tsv invalid_entropy_bc_list.tsv"
clusters="invalid_entropy_bc_list.tsv"
for f in $clusters;
do 
	bc_list_file=$res_dir/$f
	bw_dir=${bc_list_file/%.tsv/}
	mkdir -p ${bw_dir}

	N=$(wc -l < "$bc_list_file")
	for ((i=2; i<=N; i++));
	do
		bc_index=$((i -2 ))
		echo $bc_index
		bc=$(head -n $i $bc_list_file | tail -n 1)
		bc1=${bc//*[[:space:]]/}"-1"


		samtools view -h $frag_bam_file | awk -F '\t' -v OFS='\t' -v BC="$bc1" '
		{ if ($1 ~ /^@/) { print $0; next }
		if (  $1 == BC ) print $0  
		}
		'  - |  samtools sort -@ 8 -o $bw_dir/sorted_farg_bc_${bc_index}.bam 
	
		samtools index $bw_dir/sorted_farg_bc_${bc_index}.bam
		bamCoverage -b $bw_dir/sorted_farg_bc_${bc_index}.bam \
		    -o $bw_dir/sorted_farg_bc_${bc_index}.bw \
		    -bs 1 \
		    -p 8
	done

done
