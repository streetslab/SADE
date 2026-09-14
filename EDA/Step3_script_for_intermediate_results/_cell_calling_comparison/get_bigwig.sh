
source ~/python_virtuenv/atac_tobias3.13/bin/activate


pres_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'


output_dir=${pres_dir}/_cell_calling_comparison


# For each barcode 
frag_bam_dir="${output_dir}/per_bc_frag_bam2"
for f in $(ls $frag_bam_dir)
do 
	if [[ ! -f "${frag_bam_dir}/${f/%bam/bw}" ]]
	then 	
	bamCoverage -b	"${frag_bam_dir}/$f"  \
	            -o "${frag_bam_dir}/${f/%bam/bw}"  \
		    -bs 1 \
		    -p 18
	fi
done

exit 0

# For cmp_3Set
frag_bam_dir="${output_dir}/cmp_3Set"
for f in $(ls $frag_bam_dir)
do 
	if [[ ! -f "${frag_bam_dir}/${f/%bam/bw}" ]]
	then 	
	bamCoverage -b	"${frag_bam_dir}/$f"  \
	            -o "${frag_bam_dir}/${f/%bam/bw}"  \
		    -bs 1 \
		    -p 18
	fi
done

exit 0



# get bigwig reads coverage with fragments info (fragments[bed] converted to bam)
frag_bam_dir="${output_dir}/cmp_SC"
for f in $(ls $frag_bam_dir)
do 
	if [[ ! -f "${frag_bam_dir}/${f/%bam/bw}" ]]
	then 	
	samtools index "${frag_bam_dir}/$f"
	bamCoverage -b	"${frag_bam_dir}/$f"  \
	            -o "${frag_bam_dir}/${f/%bam/bw}"  \
		    -bs 1
	fi
done

frag_bam_dir2="${output_dir}/cmp_ST"
for f in $(ls $frag_bam_dir2)
do
	if [[ ! -f "${frag_bam_dir2}/${f/%bam/bw}" ]]
	then 
        samtools index "${frag_bam_dir2}/$f"
        bamCoverage -b  "${frag_bam_dir2}/$f"  \
                    -o "${frag_bam_dir2}/${f/%bam/bw}"  \
                    -bs 1
	fi
done



frag_bam_dir3="${output_dir}/cmp_DNAdebris_entropy"
for  f in $(ls $frag_bam_dir3)
do 
	if [[ ! -f "${frag_bam_dir3}/${f/%bam/bw}" ]]
	then
	samtools index "${frag_bam_dir3}/$f"
        bamCoverage -b  "${frag_bam_dir3}/$f"  \
                    -o "${frag_bam_dir3}/${f/%bam/bw}"  \
                    -bs 1
	fi
done

# normalize by CPM
frag_bam_dir3="${output_dir}/cmp_DNAdebris_entropy"
for  f in $(ls $frag_bam_dir3)
do 
	if [[ ! -f "${frag_bam_dir3}/${f/%bam/CPM.bw}" ]]
	then
        bamCoverage -b  "${frag_bam_dir3}/$f"  \
                    -o "${frag_bam_dir3}/${f/%bam/CPM.bw}"  \
		    --normalizeUsing CPM \
		    --blackListFileName /home/syyang/GitRepo/atac/ref/hg38/hg38-blacklist.bed \
		    -p 16 \
                    -bs 1
	fi
done


frag_bam_dir2="${output_dir}/cmp_ST"
for f in $(ls $frag_bam_dir2)
do

	if [[ ! -f "${frag_bam_dir3}/${f/%bam/CPM.bw}" ]]
	then
        bamCoverage -b  "${frag_bam_dir2}/$f"  \
                    -o "${frag_bam_dir2}/${f/%bam/CPM.bw}"  \
		    --normalizeUsing CPM \
		    --blackListFileName /home/syyang/GitRepo/atac/ref/hg38/hg38-blacklist.bed \
		    -p 16 \
                    -bs 1
	fi
done
