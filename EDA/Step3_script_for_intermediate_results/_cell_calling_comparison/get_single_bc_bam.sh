

bamfile='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison/filtered.bam'
pres_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'

single_bc_dir=${pres_dir}/_cell_calling_comparison/Single_bc_bam/TSS_bc_lowFRIP_bc
mkdir -p $single_bc_dir
bc="ATGTTTGGTTTAATCG"
bc="AAGCTAGGTTACTTCT"
echo $bc"-1"  > ${single_bc_dir}/barcode.txt


bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
                        -o ${single_bc_dir} \
                        -f ${single_bc_dir}/barcode.txt

mv ${single_bc_dir}/filtered.bam ${single_bc_dir}/"${bc}.bam"

 samtools index ${single_bc_dir}/"${bc}.bam"
 bamCoverage -b ${single_bc_dir}/"${bc}.bam" \
     -o ${single_bc_dir}/"${bc}.bw" \
     -bs 1 \
     -p 16
