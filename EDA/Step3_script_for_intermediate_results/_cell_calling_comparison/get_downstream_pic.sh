SCRIPT_DIR="/home/syyang/GitRepo/SADE"
source ${SCRIPT_DIR}/config.sh


res_dir="/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison"
barcode_file="${res_dir}/Union_cells_bc.txt"
frag_file="${res_dir}/union_bc_fragments.tsv.gz"

peak_dir="${res_dir}/NoDNAdebris_peaks"

    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate ${Renv_Conda} 
    # Run the R script to get the PIC matrix 
    Rscript "/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F"/scripts/count_pic_mtx.r   --barcode_file="${barcode_file}" --res_dir="${peak_dir}" --frag_file="${frag_file}"

    exit 0
peak_dir="${res_dir}/Union_cells_peaks"

    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate ${Renv_Conda} 
    # Run the R script to get the PIC matrix 
    Rscript "/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F"/scripts/count_pic_mtx.r   --barcode_file="${barcode_file}" --res_dir="${peak_dir}" --frag_file="${frag_file}"


peak_dir="${res_dir}/Entropy_cells_peaks"

    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate ${Renv_Conda} 
    # Run the R script to get the PIC matrix 
    Rscript "/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F"/scripts/count_pic_mtx.r   --barcode_file="${barcode_file}" --res_dir="${peak_dir}" --frag_file="${frag_file}"


peak_dir1="${res_dir}/TSS_cells_peaks"

    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate ${Renv_Conda} 
    # Run the R script to get the PIC matrix 
    Rscript "/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F"/scripts/count_pic_mtx.r   --barcode_file="${barcode_file}" --res_dir="${peak_dir1}" --frag_file="${frag_file}"

peak_dir2="${res_dir}/FRIP_cells_peaks"

    conda_dir=$(conda info | grep -i 'base environment' | awk '{print $4 }'  )
    source "${conda_dir}/etc/profile.d/conda.sh" 
    conda activate ${Renv_Conda} 
    # Run the R script to get the PIC matrix 
    Rscript "/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F"/scripts/count_pic_mtx.r   --barcode_file="${barcode_file}" --res_dir="${peak_dir2}" --frag_file="${frag_file}"


    conda deactivate

