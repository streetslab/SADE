#!/bin/bash

seqkit fx2tab --length --name --header-line {assembly-fasta-file} > {species}_genome_chromsize.tsv
awk -F ' '  -v OFS='\t' '{print $1, $2, $3}' {species}_genome_chromsize.tsv
# change field separator to be tab


#e.g.
mouse_genome_fa='/home/syyang/CR/refdata-cellranger-arc-GRCm39-2024-A/fasta/genome.fa'
seqkit fx2tab --length --name --header-line ${mouse_genome_fa} > mouse_genome_chromsize.tsv
awk -F ' '  -v OFS='\t' '{print $1, $2, $3}' mouse_genome_chromsize.tsv > _ 
mv _ mouse_genome_chromsize.tsv


#e.g.
mouse_genome_fa='/home/syyang/CR/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa'
seqkit fx2tab --length --name --header-line ${mouse_genome_fa} > mm10_genome_chromsize.tsv
awk -F ' '  -v OFS='\t' '{print $1, $2, $3}' mm10_genome_chromsize.tsv > _
mv _ mm10_genome_chromsize.tsv
