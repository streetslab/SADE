#!/bin/bash

sample='human_brain2_3k'
windowsize_list='3000'


#for ws in $windowsize_list; do
#	echo $ws
#
#	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/atac_fragments.tsv.gz \
#		-o ~/adipose_ln/atac/res/${sample}_WS${ws} \
#		-w ${ws} \
#		-c 'chr1' \
#		-g 'hg38'
#done


for ws in $windowsize_list; do
	echo $ws
	bash /home/syyang/GitRepo/atac/compare_cellcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws} \
		-b /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/filtered_feature_bc_matrix/barcodes.tsv 
done


#for ws in $windowsize_list; do
#	echo $ws
#	bash /home/syyang/GitRepo/atac/compare_peakcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws} \
#		 -s /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/atac_possorted_bam.bam \
#		 -g 'hg38'  
#done

