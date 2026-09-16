
awk -F '\t' -v OFS='\t' '{print $1, $3}' ~/GitRepo/SADE/ref/hg38_genome_chromsize.tsv > simplified_hg38_genome_chromsize.tsv


# sort -k1,1 -k2,2n union_bc_fragments.tsv  > union_bc_fragments.bed
bedtools bedtobam -i union_bc_fragments.bed -g simplified_hg38_genome_chromsize.tsv  > union_bc_fragments.bam 

samtools view -h union_bc_fragments.bam | awk -F '\t' -v OFS='\t'  '
NR==FNR {
    sc[$1]=$6; st[$1]=$7; dr[$1]=$10"_dnadebris"; next
  } 
  {
    if($1 ~ /^@/) { print $0; next } # Pass through header lines
    
    # Only append if the barcode exists in our metadata
    CB = "CB:Z:"$1
    SC = "SC:Z:"(sc[$1] ? sc[$1] : "NA")
    ST = "ST:Z:"(st[$1] ? st[$1] : "NA")
    DR = "DR:Z:"(dr[$1] ? dr[$1] : "0")
    
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11, CB, SC, ST, DR
       }' Union_cell_3set.tsv -   | samtools sort -@ 8 -o sorted_tagged_union_bc_fragments.bam -

samtools index sorted_tagged_union_bc_fragments.bam

