#!/bin/bash

sample='VIB_10xmultiome_2'
ws='3000'

#
	echo $ws

	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/fragments.tsv.gz \
		-o ~/adipose_ln/atac/res/${sample}_WS${ws}_test4 \
		-w ${ws} \
		-c 'chr1' \
		-g 'hg38'
#
#	echo $ws
#
#	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/fragments.tsv.gz \
#		-o ~/adipose_ln/atac/res/${sample}_WS${ws}_test5 \
#		-w ${ws} \
#		-c 'chr1 chr2 chr3 chr4 chr5' \
#		-g 'hg38'

	#echo $ws

	#bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/fragments.tsv.gz \
	#	-o ~/adipose_ln/atac/res/${sample}_WS${ws}_fullchr \
	#	-w ${ws} \
	#	-c 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22' \
	#	-g 'hg38'

#	echo $ws
#
#        bash /home/syyang/GitRepo/atac/compare_peakcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws}_fullchr \
#                 -s /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/possorted_bam.bam \
#                 -g 'hg38'  
#
#
#	bash /home/syyang/GitRepo/atac/compare_cellcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws}_fullchr \
#		-b /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv 
