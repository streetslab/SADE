
#!/bin/bash



while getopts ":o:b:c:" opt; do
  case $opt in
    o)
      comparison_dir="$OPTARG"
      comparison_dir=${comparison_dir%/}  # remove trailing slash if exists
      ;;
    b)
      base_peakset_dir="$OPTARG"
      base_peakset_dir=${base_peakset_dir%/}  # remove trailing slash if exists
      ;;
    c)
      cmp_peakset_dir="$OPTARG"
      cmp_peakset_dir=${cmp_peakset_dir%/}
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done



mkdir -p $comparison_dir

base_peak_file=$base_peakset_dir/peaks.bed
cmp_peak_file=$cmp_peakset_dir/peaks.bed



bedtools intersect -wo -a $base_peak_file -b $cmp_peak_file  >  $comparison_dir/recovered_peaks.bed  # intersection 
bedtools intersect -v -a $base_peak_file -b $cmp_peak_file  >  $comparison_dir/lost_peaks.bed  # in file-a but not in file-b
bedtools intersect -v -a $cmp_peak_file -b $base_peak_file  >  $comparison_dir/newly_discovered_peaks.bed 


