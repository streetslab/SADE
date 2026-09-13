# ATAC

CR_atac='/home/syyang/CR/cellranger-atac-2.1.0/bin/cellranger-atac'
$CR_atac count --id=VIB_10xmultiome_2 \
                        --chemistry=ARC-v1 \
                        --reference=/home/syyang/CR/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                        --fastqs=/mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/fastq/MO_V2/atac_SRR24422953 \
                        --localcores=30 \
                        --localmem=100 
