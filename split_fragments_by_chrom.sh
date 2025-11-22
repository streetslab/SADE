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

awk '{print $1}' ${fragment_file}  | uniq -c > ${chrom_counts_file}

# split fragments file by chromosome
## initialize counters
head_n=0
tail_n=0
while IFS=' ' read -r count chrom; do
    # echo " ${chrom} ,  ${count} ."
    head_n=$((head_n + count ))
    tail_n=$((count))
    
    chrom_frag_file="${chromosome_fragment_dir}/${chrom}_fragments.tsv"
    head -n ${head_n} ${fragment_file} | tail -n ${tail_n} > ${chrom_frag_file}
done < ${chrom_counts_file}


# create finish file 
finish_file=${chromosome_fragment_dir}/.split_fragments_by_chromosome.finished
touch ${finish_file}