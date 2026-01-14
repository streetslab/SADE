#!/bin/bash

filter_sh='/home/syyang/GitRepo/atac/filter_bam_with_barcodes.sh'

output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000F'
rna_subdir=${output_dir}/_rna_cell_types
atac_CBZbc_subdir=${rna_subdir}/cluster_atac_CBZbcs

results_dir=${rna_subdir}/filtered_bam

clusters='0 2 3 4 5 6 7 8 9 10'
for cluster in $clusters; do

	df_file="RNAcluster_${cluster}_df.csv"
	CBZ_file=${atac_CBZbc_subdir}/CBZ_$df_file

	result_filtered_bam_dir=${results_dir}/${df_file%_df.csv}
	#echo $result_filtered_bam_dir
	#ls $CBZ_file
	 
	mkdir -p ${result_filtered_bam_dir}

	bash ${filter_sh} -s ${output_dir}/entropy_filtered.bam  \
	                  -o ${result_filtered_bam_dir} \
	     	          -f ${CBZ_file}
done

