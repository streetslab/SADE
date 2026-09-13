#!/bin/bash

sample='human_10k_1'
windowsize_list='100 500 1000 3000 6000 10000 15000 20000 30000 50000 80000 100000 200000 400000 800000 100000 1000000'


for ws in $windowsize_list; do
	echo $ws
	if [[ $ws == 6000 ]]; then
		k=3
	else
		k=2
	fi

	if [[ ! -f ~/adipose_ln/atac/res/_${sample}_WS${ws}/entropy_cutoff.csv ]]; then
	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/10k_pbmc_ATACv1p1_nextgem_Chromium_X_fragments.tsv.gz \
		-o ~/adipose_ln/atac/res/_${sample}_WS${ws} \
		-w ${ws} \
		-c 'chr1' \
		-g 'hg38' \
		-k ${k}
	        rm ~/adipose_ln/atac/res/_${sample}_WS${ws}/fragments.tsv
		rm ~/adipose_ln/atac/res/_${sample}_WS${ws}/filtered_fragments.tsv	
	fi
done

