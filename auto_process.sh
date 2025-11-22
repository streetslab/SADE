#!/bin/bash 

_CHROMOSOME='chr1'
_WindowSize=3000 # default window size
# << GLOBAL VARIABLES >>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# environment setup
source ${SCRIPT_DIR}/config.sh
# load 
#    $PYTHON_ENV
#    $Genrich
#echo "python_path $PYTHON_ENV" # for debugging


## BEFORE ANYTHING ELSE: process options 
while  getopts "f:o:b:s:c:w:g:" opt; do
  case $opt in
    f) 
      fragment_file="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      output_dir=${output_dir%/}  # remove trailing slash if exists
      ;;
    b)
      crbarcode_file="$OPTARG"
      ;;
    s)
      bam_file="$OPTARG"
      ;;
    c)
      chromosome="$OPTARG" # has default value
      ;;
    w)
      window_size="$OPTARG"  # has default value
      ;;
    g)
      species="$OPTARG"
      # currently: 'human', 'mouse'
      ;;
    \?) 
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


# set default value if not provided
## set chromosome
if [[ -z ${chromosome} ]]; then
    chromosome=${_CHROMOSOME}
fi
## size window size 
if [[ -z ${window_size} ]]; then
    window_size=${_WindowSize}
fi
##
if [[ -z ${species} ]]; then
    echo "Please provide species with -g"
    exit 1
fi


figures_subdir="${output_dir}/figures"
mkdir -p $figures_subdir

# make sure to only execute the script on the fragment file that's within the output directory
frag_file=${output_dir}/fragments.tsv
checkpoint_file=${output_dir}/.finish_format_fragments
if [[ ! -f ${checkpoint_file} ]] ; then
    echo "Checking fragments file format for entropy calculation..."
    # sometimes fragments.tsv from CR does not have the correct record.. 
    gunzip -c $fragment_file > ${output_dir}/_fragments.tsv
    awk -F '\t' '{if (NF == 5) print $0}' ${output_dir}/_fragments.tsv > ${frag_file}
    rm ${output_dir}/_fragments.tsv 
    touch ${checkpoint_file}  # create empty file as checkpoint for fragment formatting
fi





# Record parameters used in this run
echo "Parameters used in this run:" > ${output_dir}/parameters.csv
echo "species,${species}" >> ${output_dir}/parameters.csv
echo "chromosome,${chromosome}" >> ${output_dir}/parameters.csv
echo "window_size,${window_size}" >> ${output_dir}/parameters.csv


# MODULE 0:  PREPEARATION for global variables and files to plot
## Remove # lines from the fragments file and count the number of unique fragments per barocode
checkpoint_file="${output_dir}/.finish_module_0"
if [[ ! -f ${checkpoint_file} ]] ; then
  sed '/^#/d' $frag_file| cut -f4 | sort | uniq -c > ${output_dir}/total_fragments_counts.txt
  ## Inplace trailing off the leading spaces in the output files
  sed -i 's/^[ ]*//'  ${output_dir}/total_fragments_counts.txt 
  touch ${checkpoint_file}  # create empty file as checkpoint for MODULE 0
fi




# MODULE 1:  Calculate entropies of each barcode
barcode_entropy_df_file="${output_dir}/${chromosome}_barcode_entropy_df.tsv"
if [ ! -f ${barcode_entropy_df_file} ] ; then
  bash ${SCRIPT_DIR}/calculate_entropy.sh \
          -o $output_dir \
          -g ${species} \
          -c ${chromosome}  \
          -w ${window_size}
fi 


# ++MODULE 1.5: automatically determine entropy threshold based on knee plot
auto_entropy_file="${output_dir}/entropy_cutoff.csv"
if [ ! -f ${auto_entropy_file} ] ; then
    source ${PYTHON_ENV}
    python ${SCRIPT_DIR}/find_threshold.py \
        --output_dir ${output_dir} \
        --chromosome ${chromosome}

    # Record entropy threshold into the parameters file
    awk 'NR==1' ${auto_entropy_file} >> ${output_dir}/parameters.csv
fi



# MODULE 2:  Filter fragments based on entropy thresholdded barcodes
#            and plot it overlaying with CR cell-calling labels
filtered_frag_file=${output_dir}/filtered_fragments.tsv  #check-point for MODULE 2
filtered_bc_file=${output_dir}/temp_bc_CBZ.txt  #check-point for MODULE 2
if [ ! -f ${filtered_frag_file} ] ; then 
    # retrieve entropy value from knee method 
    entropy_threshold=$(awk -F ',' 'NR==1 {print $2}'  ${auto_entropy_file}  )

    # Filter fragments based on the entropy threshold
    bash ${SCRIPT_DIR}/filter_fragments.sh \
        -o ${output_dir} \
        -c ${chromosome}
  
    # Plot [Optional]
    source ${PYTHON_ENV}
    python ${SCRIPT_DIR}/overlay_entropy_CRcelllabel_plot.py \
        --output_dir $output_dir \
        --chromosome $chromosome \
        --crbarcode_file $crbarcode_file 
fi 


# MODULE 3:  Filter BAM file with barcodes passed the entropy threshold
filtered_bam_file=${output_dir}/entropy_filtered.bam #check-point for MODULE 3
if [ ! -f ${filtered_bam_file} ] ; then
    bash ${SCRIPT_DIR}/filter_bam_with_barcodes.sh \
        -s ${bam_file} \
        -o ${output_dir} \
        -f ${filtered_bc_file}
fi


# MODULE 4:  (I) peak calling on the filtered BAM file
entropy_peak_calling_subdir="${output_dir}/peaks_entropy_filtered"
if [ ! -f ${entropy_peak_calling_subdir}/peaks.bed ] ; then
    bash ${SCRIPT_DIR}/call_peaks.sh \
    -s ${filtered_bam_file} \
    -o ${entropy_peak_calling_subdir} 
fi
# MODULE 4:  (II) peak calling on the original BAM file
peak_calling_subdir="${output_dir}/peaks"
if [ ! -f ${peak_calling_subdir}/peaks.bed ] ; then
    bash ${SCRIPT_DIR}/call_peaks.sh \
    -s ${bam_file} \
    -o ${peak_calling_subdir}
fi


# MODULE 5:  Compare the peak calling results
