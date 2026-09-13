# GEX

cellranger='/home/syyang/CR/cellranger-9.0.1/cellranger'
$cellranger count --id=VIB_10xmultiome_2_rna\
                   --chemistry=ARC-v1 \
                   --transcriptome=/home/syyang/CR/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                   --fastqs=/mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/fastq/MO_V2/rna_SRR24446480 \
                   --sample=VIB_10xmultiome_2\
                   --localcores=36 \
                   --localmem=100 \
		   --create-bam=true
