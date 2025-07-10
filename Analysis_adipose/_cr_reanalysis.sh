#!/bin/bash
sample="NK_ATAC_MAPLE013_SC"

CR_atac='/home/syyang/CR/cellranger-atac-2.1.0/bin/cellranger-atac'


output_dir="/mnt/hdd_bob/syy/adipose/atac/res/${sample}/CR_reanalysis"
cd ${output_dir}

# sort -k1,1 -k2,2n  ../filtered_fragments.tsv > filtered_fragments.sorted.tsv
# bgzip filtered_fragments.sorted.tsv
# tabix -p bed filtered_fragments.sorted.tsv.gz
#awk -v OFS='\t' '{print $1, $2, $3}' /mnt/hdd_bob/syy/adipose/atac/res/${sample}/peaks_entropy_filtered/peaks.bed > peaks_entropy_filtered.bed

${CR_atac} reanalyze --id=entropyfiltered_peaks \
                            --peaks=peaks_entropy_filtered.bed \
                            --reference=/home/syyang/CR/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                            --fragments=/mnt/hdd_bob/syy/adipose/atac/res/${sample}/CR_reanalysis/filtered_fragments_sorted.tsv.gz