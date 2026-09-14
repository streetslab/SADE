
bamfile='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison/Union_cells_peaks/filtered.bam'
pres_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'


# For cmp_SC
output_dir=${pres_dir}/_cell_calling_comparison/cmp_SC
set_list="entropy_only both none cr_only"
for set_name in $set_list
do
    awk -F '\t' -v OFS="" -v SET_NAME="$set_name"  -v prefix=CB:Z: 'NR>1 {if ($6==SET_NAME) print prefix, $1, -1}' ${pres_dir}/_cell_calling_comparison/Union_cell_3set.tsv  > ${output_dir}/${set_name}_CBZ.txt
    
    bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
                            -o ${output_dir} \
                            -f ${output_dir}/${set_name}_CBZ.txt && \
    			mv ${output_dir}/filtered.bam ${output_dir}/${set_name}_filtered.bam
done


# For cmp_ST
output_dir=${pres_dir}/_cell_calling_comparison/cmp_ST
set_list="both none tss_only"
for set_name in $set_list
do
	#set_name="entropy_only"

    awk -F '\t' -v OFS="" -v SET_NAME="$set_name"  -v prefix=CB:Z: 'NR>1 {if ($7==SET_NAME) print prefix, $1, -1}' ${pres_dir}/_cell_calling_comparison/Union_cell_3set.tsv  > ${output_dir}/${set_name}_CBZ.txt
    # ^ keep a column for barocde with '-1' appendix for filtering fragment file
    
    bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
                            -o ${output_dir} \
                            -f ${output_dir}/${set_name}_CBZ.txt && \
    			mv ${output_dir}/filtered.bam ${output_dir}/${set_name}_filtered.bam
done




# For cmp_3Set
output_dir=${pres_dir}/_cell_calling_comparison/cmp_3Set
mkdir -p $output_dir
set_list=("common" "CR unique" "SADE unique" "TSS unique" "SADE+CR, not TSS" "SADE+TSS, not CR" "CR+TSS, not SADE")
#for set_name in $set_list
for i in {1..6}
do
	set_name=${set_list[$i]}
	echo $set_name

    awk -F '\t' -v OFS="" -v SET_NAME="$set_name"  -v prefix=CB:Z: 'NR>1 {if ($2==SET_NAME) print prefix, $1 }' ${pres_dir}/_cell_calling_comparison/_union_cell_final_info.csv  > ${output_dir}/${set_name//[" ",]/_}_CBZ.txt
    
    bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
                            -o ${output_dir} \
                            -f ${output_dir}/${set_name//[" ",]/_}_CBZ.txt && \
    			mv ${output_dir}/filtered.bam ${output_dir}/${set_name//[" ",]/_}_filtered.bam
done

