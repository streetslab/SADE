#!/bin/bash

gzip -cd  cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz  > /mnt/files/syy/adipose/multiom/Arionas_data_May2023/Analysis/SCVI_integration/rna_barcodes.tsv
gzip -cd  cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz  > /mnt/files/syy/adipose/multiom/Arionas_data_May2023/Analysis/SCVI_integration/atac_barcodes.tsv

paste atac_barcodes.tsv  rna_barcodes.tsv | awk -v OFS='\t' '{print , }' >  atac_rna_barcodes_map.tsv
