
output_dir='/home/syyang/adipose_ln/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison'
file=$output_dir/Union_cell_3set_with_RNA_cluster_and_DNAdebrisflag.tsv

awk -F '\t' -v OFS='\t'  '
NR==1 {print $0, "_3set_identified_by_atac"; next }
{
 if ($3=="True" && $4=="True" && $5=="True") {
	print $0, "common"}
 else if ($3=="False" && $4=="True" && $5=="True") {
	 print $0, "SADE+TSS, not CR"}
 else if ($3=="True" && $4=="False" && $5=="True") {
	 print $0, "CR+TSS, not SADE"}
 else if ($3=="True" && $4=="True" && $5=="False") {
	 print $0, "SADE+CR, not TSS"}
 else if ($3=="False" && $4=="False" && $5=="True") {
	       print $0, "TSS unique"}
 else if ($3=="False" && $4=="True" && $5=="False") {
	       print $0, "SADE unique"}
 else if ($3=="True" && $4=="False" && $5=="False") {
	       print $0, "CR unique"}
 }' $file >  $output_dir/Union_cell_3set_all_info.tsv
