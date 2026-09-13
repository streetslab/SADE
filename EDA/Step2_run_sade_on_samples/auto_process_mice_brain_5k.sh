#!/bin/bash

sample='mice_brain_5k'
#windowsize_list='200 300 350 400 450 500 100 700 1000 2000 3000 5000 7000'
#windowsize_list='1500 1700 1900 2300 2500 2800 3100 3300 3500 3700 3900 40000 10000 30000 50000 80000 100000'
#windowsize_list='500 1500 5000 10000 30000 100000 1000000'
windowsize_list='3000'


source  ~/python_virtuenv/mrvi_3.11/bin/activate

for ws in $windowsize_list; do
	echo $ws

	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/fragments.tsv.gz \
		-o ~/adipose_ln/atac/res/${sample}_WS${ws} \
		-c 'chr1' \
		-w ${ws} \
		-g 'mm10'
	bash /home/syyang/GitRepo/atac/compare_cellcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws} \
		-b /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/filtered_feature_bc_matrix/barcodes.tsv \
	

	bash /home/syyang/GitRepo/atac/compare_peakcalling.sh -d ~/adipose_ln/atac/res/${sample}_WS${ws} \
		-s /mnt/hdd_bob/syy/adipose/atac/cr_goldstand/${sample}/possorted_bam.bam \
		-g 'mm10'


done






