# atac
Use entropy to do cell-calling on scATAC-seq data. 
Benefit: Increase peak signal detection sensitivity.  

## 
auto_process.sh calculate per barcode entropy and auto-threshold entropy


## 
compare_cellcalling.sh compares CellRaner's (post-peak) cell-calling method v.s. Entropy (pre-peak) cell-calling


## 
compare_peakcalling.sh compares peak results before- v.s. post- Entropy filtering


##
downstream_analysis.sh compares peak featuers (before- v.s. post- Entropy filtering)'s ability for cell-type discovery 




### STEP 1. 
bash auto_process.sh 
### STEP 2. 
bash compare_cellcalling.sh
### STEP 3. 
bash compare_peakcalling.sh
### STEP 4. 
bash downstream_analysis.sh
