#!/bin/bash

sample='mice_brain_5k'
windowsize_list='100 500 1000 3000 6000 10000 15000 20000 30000 50000 80000 100000 200000 400000 800000 100000 1000000'


source  ~/python_virtuenv/mrvi_3.11/bin/activate

for ws in $windowsize_list; do
	echo $ws
        if [[ $ws == 500 || $ws == 1000 ]]; then
                k=3
        else
                k=2
        fi

	echo $k

	if [[ ! -f ~/adipose_ln/atac/res/_${sample}_WS${ws}/entropy_cutoff.csv ]]; then
		bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/fragments.tsv.gz \
		-o ~/adipose_ln/atac/res/_${sample}_WS${ws} \
		-c 'chr1' \
		-w ${ws} \
		-g 'mm10' \
		-k ${k}
		rm ~/adipose_ln/atac/res/_${sample}_WS${ws}/fragments.tsv
		rm ~/adipose_ln/atac/res/_${sample}_WS${ws}/filtered_fragments.tsv
	fi
done






