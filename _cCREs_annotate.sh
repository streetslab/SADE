#!/bin/bash

output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000F'
output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000_test_norm'




cCREs_annotation_file='/home/syyang/GitRepo/atac/ref/hg38/GRCh38-cCREs.bed'

comparison_subdir=${output_dir}/peaks_comparison

newly_discover_peaks=$comparison_subdir/newly_discovered_peaks.bed 
lost_peaks=$comparison_subdir/lost_peaks.bed 


cCRE_subdir=${output_dir}/cCRE_annotation
mkdir -p ${cCRE_subdir}
bedtools intersect -wao -a $newly_discover_peaks -b $cCREs_annotation_file > ${cCRE_subdir}/cCRE_annot_newly_discovered_peaks.bed  #intersection
bedtools intersect -wao -a $lost_peaks -b $cCREs_annotation_file > ${cCRE_subdir}/cCRE_annot_lost_peaks.bed  
