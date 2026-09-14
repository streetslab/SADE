#!/bin/bash

# Should be executed within the output_dir from running sade.sh
# ie. within the VIB_10xmultiome_2_WS3000F dir 

# Uncomment to run
cp ../fragments.tsv  . # from the original fragments.tsv file 
bgzip fragments.tsv
tabix --preset bed fragments.tsv.gz
#
### Run the R-script  get_tss_w_archr.r
# Activate conda env that has R
Renv_Conda='seuratv5'
conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
source "${conda_dir}/etc/profile.d/conda.sh" 
conda activate ${Renv_Conda} 

Rscript get_tss_w_archr.r

head -n 1 QualityControl/VIB_10xmultiome_2_Metadata.tsv  > TSS_passed_bc.tsv  # write ArchR header
awk -F '\t' -v OFS='\t' '{if ($6>4) print $0}'   QualityControl/VIB_10xmultiome_2_Metadata.tsv | sort -r -k6,6 -g  >> TSS_passed_bc.tsv
# cat  QualityControl/VIB_10xmultiome_2_Metadata.tsv  | sort -r -k6,6 -g  >> TSS_passed_bc.tsv

awk -F'\t' 'NR>1 {print $1}'  TSS_passed_bc.tsv > TSS_pass_bc.csv
sed -i 's/"//g' TSS_pass_bc.csv 
sed -i 's/VIB_10xmultiome_2#//g' TSS_pass_bc.csv 
awk -v OFS='' -v prefix='CB:Z:' '{print prefix, $1}' TSS_pass_bc.csv > TSS_pass_bc_CBZ.txt



rm fragments.tsv.gz