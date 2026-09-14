

output_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'
script_dir=$output_dir/scripts
res_dir=$output_dir/_cell_calling_comparison

# Each nucleus calling method's nuclei set's reads peaks VS raw peaks
bash ${script_dir}/cmp_peaksets.sh -o $res_dir/FRIPpeaks_vs_Rawpeaks \
	-b  $res_dir/Rawpeaks \
	-c  $res_dir/FRIP_cells_peaks


bash ${script_dir}/cmp_peaksets.sh -o $res_dir/TSSpeaks_vs_Rawpeaks \
	-b $res_dir/Rawpeaks \
	-c  $res_dir/TSS_cells_peaks

bash ${script_dir}/cmp_peaksets.sh -o $res_dir/Entropypeaks_vs_Rawpeaks \
	-b $res_dir/Rawpeaks \
	-c $res_dir/Entropy_cells_peaks

exit 0


# Common vs Raw
bash ${script_dir}/cmp_peaksets.sh -o $res_dir/Commonpeaks_vs_Rawpeaks \
	-b  $res_dir/Rawpeaks \
	-c  $res_dir/Commoncells_peaks


#  VS common nuclei's inferred peaks
bash ${script_dir}/cmp_peaksets.sh -o $res_dir/FRIPpeaks_vs_Commonpeaks \
	-b $res_dir/Commoncells_peaks \
	-c $res_dir/FRIP_cells_peaks

bash ${script_dir}/cmp_peaksets.sh -o $res_dir/TSSpeaks_vs_Commonpeaks \
	-b $res_dir/Commoncells_peaks \
	-c $res_dir/TSS_cells_peaks

bash ${script_dir}/cmp_peaksets.sh -o $res_dir/Entropypeaks_vs_Commonpeaks \
	-b $res_dir/Commoncells_peaks \
	-c $res_dir/Entropy_cells_peaks


# Entropy peaks vs TSS peaks
bash ${script_dir}/cmp_peaksets.sh -o $res_dir/Entropypeaks_vs_TSSpeaks \
	-b $res_dir/TSS_cells_peaks \
	-c $res_dir/Entropy_cells_peaks




