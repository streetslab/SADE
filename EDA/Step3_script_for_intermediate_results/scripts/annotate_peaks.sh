
SCRIPT_DIR="/home/syyang/GitRepo/atac"
source ${SCRIPT_DIR}/config.sh


while getopts ":d:g:" opt; do
  case $opt in
    d)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    g)
      genome_name="$OPTARG"
      ;;
  esac
done



peaks=$output_dir/peaks.bed
newly_discovered_peaks=$output_dir/newly_discovered_peaks.bed 
lost_peaks=$output_dir/lost_peaks.bed


annotation_subdir="${output_dir}/annotation"
mkdir -p $annotation_subdir
annotated_peaks=$annotation_subdir/annotated_peaks.bed
annotated_newly_discovered_peaks=$annotation_subdir/annotated_newly_discovered_peaks.bed
annotated_lost_peaks=$annotation_subdir/annotated_lost_peaks.bed


# annotate peaks differences 
if [[ -f $peaks ]]
then 
	annotatePeaks.pl $peaks $genome_name -annStats  $annotation_subdir/"annStats_peaks.txt" > $annotated_peaks
fi
if [[ -f $newly_discovered_peaks ]]
then 
	annotatePeaks.pl $newly_discovered_peaks $genome_name -annStats  $annotation_subdir/"annStats_newly_discovered_peaks.txt" > $annotated_newly_discovered_peaks
fi
if [[ -f $lost_peaks ]]
then
	annotatePeaks.pl $lost_peaks $genome_name -annStats  $annotation_subdir/"annStats_lost_peaks.txt" > $annotated_lost_peaks
fi 

