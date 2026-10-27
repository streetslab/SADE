#!/bin/bash

Msg="Usage: $0 -o output_dir"

while getopts "o:" opt; do
  case ${opt} in
    o)
      output_dir=$OPTARG
      output_dir=${output_dir%/} # Remove trailing slash if present
      ;;
    *)
      echo $Msg  && exit 1
      ;;
  esac
done

if [ $# -eq 0 ] ; then
    echo $Msg && exit 1
fi

fragment_file=${output_dir}/fragments.tsv

chromosome_fragment_dir=${output_dir}/chromosome_fragments

if [ ! -d ${chromosome_fragment_dir} ] ; then
    mkdir -p ${chromosome_fragment_dir}
fi



# calculate rows for each chromosome 
chrom_counts_file="${chromosome_fragment_dir}/chromsome_row_counts.txt"

#awk '{print $1}' ${fragment_file}  | uniq -c > ${chrom_counts_file}
cut -f1 ${fragment_file}  | uniq -c > ${chrom_counts_file}

# split fragments file by chromosome
# ## initialize counters
# head_n=0
# tail_n=0
# while IFS=' ' read -r count chrom; do
#     # echo " ${chrom} ,  ${count} ."
#     chrom_frag_file="${chromosome_fragment_dir}/${chrom}_fragments.tsv"

#     head_n=$((head_n + count ))
#     tail_n=$((count))
#     # Approach1: using head and tail
#     head -n ${head_n} ${fragment_file} | tail -n ${tail_n} > ${chrom_frag_file}
#     # Approach2: Use sed to extract -- slow [Surprisingly sed is slower than head/tail combination]
#     # sed -n "$((head_n - tail_n + 1)),$((head_n))p" ${fragment_file} > ${chrom_frag_file}

# done < ${chrom_counts_file}

# Alternative approach using awk >>> Gemini 3 code
awk -v out_dir="${chromosome_fragment_dir}" '
    # 1. Read the counts file first (NR==FNR)
    NR==FNR {
        chrom_counts[NR] = $1  # Store the count
        chrom_names[NR] = $2   # Store the name
        total_chroms = NR
        next
    }
    # 2. Process the main fragment file
    FNR == 1 { 
        current_idx = 1
        lines_written = 0
        # Set first output filename
        current_file = out_dir "/" chrom_names[current_idx] "_fragments.tsv"
    }

    {
        # Write the current line to the current file
        print > current_file
        lines_written++
        if (lines_written >= chrom_counts[current_idx]) {
            close(current_file) 
            
            # Setup for next chunk
            current_idx++
            lines_written = 0
            
            # Stop if we run out of chromosomes defined in counts file
            if (current_idx > total_chroms) exit
            
            current_file = out_dir "/" chrom_names[current_idx] "_fragments.tsv"
        }
    }
' "${chrom_counts_file}" "${fragment_file}"
# Alternative approach using awk <<< Gemini 3 code


# create finish file 
finish_file=${chromosome_fragment_dir}/.split_fragments_by_chromosome.finished
touch ${finish_file}