#!/bin/bash

#sample='BJ_Cell_LAND'
#sample='BJ_CELL_xSDS'
#sample='BJ_Cell_ATAC'
windowsize_list='3000'

sample_list="BJ_Cell_ATAC BJ_CELL_xSDS BJ_Cell_LAND"
#sample_list="BJ_Cell_LAND"

for sample in $sample_list; do 

for ws in $windowsize_list; do
	echo $ws

	bash /home/syyang/GitRepo/atac/auto_process.sh -f /mnt/hdd_bob/syy/adipose/atac/defnd_seq/${sample}/fragments.tsv.gz \
		-o ~/adipose_ln/atac/res/${sample}_F \
		-w ${ws} \
		-g 'hg38'
done


done
