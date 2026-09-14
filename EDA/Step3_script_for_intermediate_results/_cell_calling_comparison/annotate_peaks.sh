
res_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison'


bash ../scripts/annotate_peaks.sh -d $res_dir/Commonpeaks_vs_Rawpeaks -g hg38
bash ../scripts/annotate_peaks.sh -d $res_dir/Commoncells_peaks -g hg38

bash ../scripts/annotate_peaks.sh -d $res_dir/Entropypeaks_vs_TSSpeaks -g hg38

bash ../scripts/annotate_peaks.sh -d $res_dir/Rawpeaks -g hg38 

bash ../scripts/annotate_peaks.sh -d $res_dir/Entropy_cells_peaks -g hg38 
bash ../scripts/annotate_peaks.sh -d $res_dir/Entropypeaks_vs_Rawpeaks -g hg38 
bash ../scripts/annotate_peaks.sh -d $res_dir/FRIP_cells_peaks -g hg38 
bash ../scripts/annotate_peaks.sh -d $res_dir/TSS_cells_peaks -g hg38 
bash ../scripts/annotate_peaks.sh -d $res_dir/FRIPpeaks_vs_Rawpeaks  -g hg38 
bash ../scripts/annotate_peaks.sh -d $res_dir/TSSpeaks_vs_Rawpeaks -g hg38 
