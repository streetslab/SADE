source ~/python_virtuenv/atac_tobias3.13/bin/activate


output_dir="/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison"

frag_bam_file=$output_dir/sorted_tagged_union_bc_fragments.bam

pres_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F'
bw_dir=${pres_dir}/_cell_calling_comparison/per_bc_frag_bam2

mkdir -p $bw_dir

cd $bw_dir

awk -F '\t' 'NR>1 {print $1}' ${pres_dir}/_cell_calling_comparison/Union_cell_3set.tsv  | while read -r bc1
do
	# header 
	samtools view -H $frag_bam_file >  "${bc1}_fragments.sam"
done

# Split the fragments into separate files
samtools view   $frag_bam_file |  awk  -F '\t' '{
	filename= $1"_fragments.sam";
	print $0 >> filename;
	close(filename);
}' - 

awk -F '\t' 'NR>1 {print $1}' ${pres_dir}/_cell_calling_comparison/Union_cell_3set.tsv  | while read -r bc1
do
	samtools view -bS "${bc1}_fragments.sam" > "${bc1}_fragments.bam"
	samtools index "${bc1}_fragments.bam"
done





#awk -F '\t' 'NR>1 {print $1}' ${pres_dir}/_cell_calling_comparison/Union_cell_3set.tsv  | while read -r bc1
#do
#	((count++))
#	if [[ count -gt 98 ]]  #Remove
#	then                   #Remove
#	samtools view -h $frag_bam_file | awk -F '\t' -v OFS='\t' -v BC="$bc1" '
#	{ if ($1 ~ /^@/) { print $0; next }
#	if (  $1 == BC ) print $0  
#	}
#	'  - |  samtools sort -@ 16 -o $bw_dir/${bc1}_fragments.bam 
#	fi                     #Remove
#done
#
