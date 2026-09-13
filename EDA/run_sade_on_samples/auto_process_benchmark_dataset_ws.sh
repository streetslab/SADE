#!/bin/bash

sample='VIB_10xmultiome_2'
windowsize_list='100 500 1000 3000 6000 10000 15000 20000 30000 50000 80000 100000 200000 400000 800000 100000 1000000'



for ws in $windowsize_list; do
	echo $ws

	if [[ ! -f ~/adipose_ln/atac/res/_${sample}_WS${ws}/entropy_cutoff.csv ]]; then

	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/fragments.tsv.gz \
		-o ~/adipose_ln/atac/res/_${sample}_WS${ws} \
		-w ${ws} \
		-c 'chr1' \
		-g 'hg38'
	rm ~/adipose_ln/atac/res/_${sample}_WS${ws}/fragments.tsv
	rm ~/adipose_ln/atac/res/_${sample}_WS${ws}/filtered_fragments.tsv
	fi 
done





