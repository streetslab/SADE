#!/bin/bash
this_github_analysis_dir='/home/syyang/GitRepo/atac/EDA'
output_dir='/mnt/hdd_bob/syy/adipose/atac/res/test_run'

# Prepare script
cp -r  "${this_github_analysis_dir}/Step3_script_for_intermediate_results/scripts" "${output_dir}/".
cp -r "${this_github_analysis_dir}/Step3_script_for_intermediate_results/_ArchR_TSS" "${output_dir}/".
cp -r  "${this_github_analysis_dir}/Step3_script_for_intermediate_results/_cell_calling_comparison" "${output_dir}/"