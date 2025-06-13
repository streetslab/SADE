bash process.sh -f /mnt/hdd_bob/syy/adipose/atac/cr_atac_results/NK_ATAC_MAPLE013_CS/outs/fragments.tsv.gz   -o res/NK_ATAC_MAPLE013_CS/ -b /mnt/hdd_bob/syy/adipose/atac/cr_atac_results/NK_ATAC_MAPLE013_CS/outs/filtered_peak_bc_matrix/barcodes.tsv -s /mnt/hdd_bob/syy/adipose/atac/cr_atac_results/NK_ATAC_MAPLE013_CS/outs/possorted_bam.bam  > res/log_NK_ATAC_MAPLE013_CS 



 bash process.sh -f /mnt/hdd_bob/syy/adipose/atac/cr_atac_results/AB_ATAC_FL_TAM_SQ/outs/fragments.tsv.gz  -o res/AB_ATAC_FL_TAM_SQ/ -b /mnt/hdd_bob/syy/adipose/atac/cr_atac_results/AB_ATAC_FL_TAM_SQ/outs/filtered_peak_bc_matrix/barcodes.tsv  -s /mnt/hdd_bob/syy/adipose/atac/cr_atac_results/AB_ATAC_FL_TAM_SQ/outs/possorted_bam.bam
 

 ## run compare
 bash ~/GitRepo/atac/compare.sh -d /mnt/hdd_bob/syy/adipose/atac/res/AB_ATAC_FL_TAM_SQ/ -g hg38


 bash GitRepo/atac/compare.sh   -d /mnt/hdd_bob/syy/adipose/atac/res/pbmc_1k/ -g hg38
