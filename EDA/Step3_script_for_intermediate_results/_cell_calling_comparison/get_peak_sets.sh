
SCRIPT_DIR='/home/syyang/GitRepo/SADE'
source ${SCRIPT_DIR}/config.sh

bamfile='/mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/VIB_10xmultiome_2/outs/possorted_bam.bam'
pres_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'
raw_peaks=${pres_dir}/peaks
entropy_peaks=${pres_dir}/peaks_entropy_filtered


output_dir=${pres_dir}/_cell_calling_comparison



Raw_peak_dir=${output_dir}/Rawpeaks
Entropy_peak_dir=${output_dir}/Entropy_cells_peaks
FRIP_peak_dir=${output_dir}/FRIP_cells_peaks
TSS_peak_dir=${output_dir}/TSS_cells_peaks


# mkdir -p ${FRIP_peak_dir} ${TSS_peak_dir} 
# mkdir -p ${Entropy_peak_dir} ${Raw_peak_dir}
# cp /mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/Entropy_filtered_bc_CBZ.txt  ${Entropy_peak_dir}/Entropy_filtered_bc_CBZ.txt
# cp ${bamfile} ${Raw_peak_dir}/possorted_bam.bam

# # for Raw peaks that uses all reads 
# bash ${SCRIPT_DIR}/call_peaks.sh -s ${Raw_peak_dir}/possorted_bam.bam \
# 	-o ${Raw_peak_dir}
# bedtools subtract -a ${Raw_peak_dir}/peaks_w_blacklistregion.bed  -b ${SCRIPT_DIR}/ref/hg38/hg38-blacklist.bed  > ${Raw_peak_dir}/peaks.bed


# # for FRIP filtered cells' peaks
# cp '/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_CR_FRIP/barcodes.tsv' ${FRIP_peak_dir}/CR_barcodes.tsv
# awk -v OFS='' -v prefix='CB:Z:' '{print prefix, $1}'  ${FRIP_peak_dir}/CR_barcodes.tsv > ${FRIP_peak_dir}/CR_filtered_bc_CBZ.txt

# bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
# 	                -o ${FRIP_peak_dir}  \
# 			-f ${FRIP_peak_dir}/CR_filtered_bc_CBZ.txt


# bash ${SCRIPT_DIR}/call_peaks.sh -s ${FRIP_peak_dir}/filtered.bam \
# 	-o ${FRIP_peak_dir}
# bedtools subtract -a ${FRIP_peak_dir}/peaks_w_blacklistregion.bed  -b ${SCRIPT_DIR}/ref/hg38/hg38-blacklist.bed  > ${FRIP_peak_dir}/peaks.bed


# For TSS filtered cells' peaks
cp ${pres_dir}/_ArchR_TSS/TSS_pass_bc_CBZ.txt  ${TSS_peak_dir}/TSS_filtered_bc_CBZ.txt

bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
	                -o ${TSS_peak_dir} \
			-f ${TSS_peak_dir}/TSS_filtered_bc_CBZ.txt

bash /home/syyang/GitRepo/atac/call_peaks.sh -s ${TSS_peak_dir}/filtered.bam \
	-o  ${TSS_peak_dir} 
bedtools subtract -a ${TSS_peak_dir}/peaks_w_blacklistregion.bed  -b ${SCRIPT_DIR}/ref/hg38/hg38-blacklist.bed  > ${TSS_peak_dir}/peaks.bed

# # for Entropy filtered cells' peaks 
# cp ${pres_dir}/Entropy_filtered_bc_CBZ.txt  ${Entropy_peak_dir}/Entropy_filtered_bc_CBZ.txt
# bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
# 	                -o ${Entropy_peak_dir}  \
# 			-f ${Entropy_peak_dir}/Entropy_filtered_bc_CBZ.txt

# bash ${SCRIPT_DIR}/call_peaks.sh -s ${Entropy_peak_dir}/filtered.bam \
# 	-o ${Entropy_peak_dir}
# bedtools subtract -a ${Entropy_peak_dir}/peaks_w_blacklistregion.bed  -b ${SCRIPT_DIR}/ref/hg38/hg38-blacklist.bed  > ${Entropy_peak_dir}/peaks.bed


# --- Infer peaks with union cells 
Common_peak_dir=${output_dir}/Union_cells_peakss
mkdir -p ${Common_peak_dir}
awk -F '\t' -v OFS='\t' '{if ($8 != "none") print $0}'  ${output_dir}/Union_cell_3set_all_info.tsv >  ${output_dir}/Union_cell_3set.tsv
awk -F '\t' -v OFS="" -v prefix=CB:Z: 'NR>1 {print prefix, $1, -1}' ${output_dir}/Union_cell_3set.tsv  > ${output_dir}/Union_cells_bc_CBZ.txt
cp ${output_dir}/Union_cells_bc_CBZ.txt   ${Common_peak_dir}/Union_cells_bc_CBZ.txt

bash ${pres_dir}/scripts/filter_bam_w_bc.sh -s ${bamfile} \
	                -o ${Common_peak_dir}  \
			-f ${Common_peak_dir}/Union_cells_bc_CBZ.txt

bash ${SCRIPT_DIR}/call_peaks.sh -s ${Common_peak_dir}/filtered.bam \
	-o ${Common_peak_dir}
bedtools subtract -a ${Common_peak_dir}/peaks_w_blacklistregion.bed  -b ${SCRIPT_DIR}/ref/hg38/hg38-blacklist.bed  > ${Common_peak_dir}/peaks.bed

#
#
# Remove those filtered bam files to save space
rm -rf ${FRIP_peak_dir}/filtered.bam 
rm -f ${TSS_peak_dir}/filtered.bam 
rm -f ${Entropy_peak_dir}/filtered.bam 
# rm -f ${Common_peak_dir}/filtered.bam