#!/bin/bash
this_github_analysis_dir='/home/syyang/GitRepo/SADE/EDA'
output_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'

# Prepare script
cp -r  "${this_github_analysis_dir}/Step3_script_for_intermediate_results/scripts" "${output_dir}/".
cp -r "${this_github_analysis_dir}/Step3_script_for_intermediate_results/_ArchR_TSS" "${output_dir}/".
cp -r  "${this_github_analysis_dir}/Step3_script_for_intermediate_results/_cell_calling_comparison" "${output_dir}/"
cp -r "${this_github_analysis_dir}/Step3_script_for_intermediate_results/_CR_FRIP" "${output_dir}/"