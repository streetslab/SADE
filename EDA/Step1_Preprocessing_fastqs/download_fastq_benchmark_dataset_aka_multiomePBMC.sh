#!/bin/bash

# Download fastqs for the multiomePBMC dataset used in the manuscript.
fastq-dump --split-files --gzip SRR24422953 # ATAC fastqs for MO_V2 (aka 'multiomePBMC', aka 'VIB_10xmultiome_2' for the sample folder name) used in the manuscript. 
fastq-dump --split-files --gzip SRR24446480 # RNA fastqs for MO_V2 (aka 'multiomePBMC', aka 'VIB_10xmultiome_2' for the sample folder name) used in the manuscript. 

## CellRanger was run on the fastqs to generate the output files used in the manuscript.
## Resulting fragments.tsv was used for SADE

# Download the CellRanger output files for the other 4 datasets from 10x 
# 1. Mouse brain 5k data.
# from: 
# https://www.10xgenomics.com/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-1-0-0
# 2. Mouse brain cryo-preserved 5k data.
# from:
# https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_E18_brain_cryo_5k/atac_v1_E18_brain_cryo_5k_web_summary.html
# 3.Human PBMC 10k data.
# from:
# https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v1-1-chromium-x-1-1-standard
# 4. Human brain 3k data.
# from:
# https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-1-0-0
## For each dataset, fragments.tsv was used for SADE