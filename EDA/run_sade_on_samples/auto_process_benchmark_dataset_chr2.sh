#!/bin/bash

sample='VIB_10xmultiome_2'
windowsize_list='500 1000 2000 3000'
windowsize_list='3000'



for ws in $windowsize_list; do
	echo $ws

	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/fragments.tsv.gz \
		-o ~/adipose_ln/atac/res/${sample}_WS${ws}_chr2 \
		-w ${ws} \
		-c 'chr2' \
		-g 'hg38'
done

#for ws in $windowsize_list; do
#	# for record 
#	grep "depth"  /home/syyang/GitRepo/atac/calculate_entropy.py  > ~/adipose_ln/atac/res/${sample}_WS${ws}_chr2/calculate_entrop.log
#done
#
#
#
#for ws in $windowsize_list; do
#	echo $ws
#	bash /home/syyang/GitRepo/atac/compare_cellcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws}_chr2 \
#		-b /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv 
#done
#
#
for ws in $windowsize_list; do
	echo $ws
	bash /home/syyang/GitRepo/atac/compare_peakcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws}_chr2 \
		 -s /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/possorted_bam.bam \
		 -g 'hg38'  
done



for ws in $windowsize_list; do
	echo $ws

	bash /home/syyang/GitRepo/atac/downstream_analysis.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws}_chr2  \
		-b ~/adipose_ln/atac/res/${sample}_WS${ws}_chr2/bc_pass_entropy.tsv

done

