# SADE
**S**hannon-entropy to filter noisy reads on single-cell **A**TAC-seq for sensitive accessible region **de**tection

------
## Run module-1  [auto_process.sh]
> [!note]
**auto_process.sh** : Calculate per barcode entropy, auto-threshold entropy and filter fragments based on barcodes.

### Calculate per barcode entropy, auto-threshold entropy and filter fragments based on barcodes.

- **System Requirement**

[ripgrep](https://github.com/BurntSushi/ripgrep)


- **Environment setup** to run **Module-1** **auto_process.sh**

> [!Required]
A python3 virtual environemnt is all you need for this module to calculate per-barcode entropy and autothresholding to filter good quality barcodes. 
1. Create a python virtual environement with required libraries
```
python -m venv you_venv1_name
```
2. Install below libraries within the virtual environment
```
# Activate your python virtual environment
source absolute_path_to_you_venv1_name/bin/activate

# install required python packages within the virtual environment
pip install -r path_to_this_dir/EnvironmentSetup/requirements_auto_process.txt
```
3. Locate your virtual environement and modify _config.sh_ by adding a line below
```
echo "export PYTHON_ENV='absolute_path_to_you_venv1_name/bin/activate'" > config.sh
```

- **Usage** of **Module-1** command:
```
sample='VIB_10xmultiome_2'
ws=3000
bash path_to_this_dir/auto_process.sh -o your_desired_output_directory \
      -f /mnt/hdd_bob/syy/adipose/atac/protocol_benchmark/cr_results/atac/${sample}/outs/fragments.tsv.gz \
      -w ${ws} \
      -c 'chr1' \
      -g 'hg38'
  ```   
  
- Command **Explanation**:
  * **Input**:   
          **-o <output_dir. Required>**
          **-f <fragments_file. Required>**  
          **-g <genome_used_for_read_mapping_that_resulted_provided_fragments_file.  Required: 'hg38', 'mm10' etc.>**  
          [-c <chromosome>. Chromosome used to calculate entropy. Default: largest chromosome, chr1 ]  
          [-w <window_size>. Windowsize on genome to look for Tn5 insertion frequencies for entropy calculation. Default: 3000]  

  * **Output**:   
          **fragments.tsv**  A copy of input frament file.   
          **filtered_fragments.tsv**  Fragments corresponding to quality nuclei by entropy criterion. This file should be used for downstream analysis in replace of fragments file.   
          **entropy_filtered_bc_df.tsv**  File with filtered barcodes as row indices, columns having entropy calculation metrics. First row has column names.
          **figures (subfolder)**
    
------ 
## Run module-2 [compare_cellcalling.sh] TODO
### Compares CellRaner's (post-peak) cell-calling method v.s. Entropy (pre-peak) cell-calling

------
## Run module-3 [compare_peakcalling.sh] TODO
### Compares CellRaner's (post-peak) cell-calling method v.s. Entropy (pre-peak) cell-calling

------
## Run module-4 [downstream_analysis.sh] 
### Compares peak featuers (before- v.s. post- Entropy filtering)'s ability for cell-type discovery 
