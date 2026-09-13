#!/bin/bash

sample='mouse_cryo_brain_5k'
windowsize_list='3000'


source  ~/python_virtuenv/mrvi_3.11/bin/activate

for ws in $windowsize_list; do
	echo $ws

	#bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/atac_v1_E18_brain_cryo_5k_fragments.tsv.gz \
	#	-o ~/adipose_ln/atac/res/${sample}_WS${ws} \
	#	-c 'chr1' \
	#	-w ${ws} \
	#	-g 'mm10'
	
	bash /home/syyang/GitRepo/atac/compare_cellcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws} \
		-b /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/filtered_peak_bc_matrix/barcodes.tsv \
	

#	bash /home/syyang/GitRepo/atac/compare_peakcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws} \
#		-s /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/possorted_bam.bam \
#		-g 'mm10'


done






