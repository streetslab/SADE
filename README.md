# atac
Use entropy to do cell-calling on scATAC-seq data. 
Benefit: Increase peak signal detection sensitivity.  


## [required]
auto_process.sh Calculate per barcode entropy and auto-threshold entropy  
> [!NOTE]
> Environment setup to run this module check [module-1 environment setup](https://github.com/Irisapo/atac/tree/c8691402345fdf6221e42f532435b6ca52d85294/EnvironmentSetup)

## 
compare_cellcalling.sh Compares CellRaner's (post-peak) cell-calling method v.s. Entropy (pre-peak) cell-calling


## 
compare_peakcalling.sh Compares peak results before- v.s. post- Entropy filtering


##
downstream_analysis.sh Compares peak featuers (before- v.s. post- Entropy filtering)'s ability for cell-type discovery 




### STEP 1. 
bash auto_process.sh 
### STEP 2. 
bash compare_cellcalling.sh
### STEP 3. 
bash compare_peakcalling.sh
### STEP 4. 
bash downstream_analysis.sh
