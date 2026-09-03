# SADE
**S**h**a**nnon-entropy **de**tection of quality nuclei in scATAC-seq.

**SADE** quantifies the complexity of genome-wide chromatin accessibility for each droplet(-barcode) with entropy, and auto-threshold entropy values to identify droplets with quality nuclei, while filtering out empty droplets, droplets with damaged nuclei, as well as droplets with cellular-debris.  


## System Requirements
[ripgrep](https://github.com/BurntSushi/ripgrep)


## Environment Setup

To run **auto_process.sh**, a Python 3 virtual environment is required to calculate per-(droplet-)barcode entropy and perform auto-thresholding.
1. Create a python virtual environement with required libraries
```bash
python -m venv you_venv1_name
```
2. Activate your python virtual environment
```bash \
source absolute_path_to_you_venv1_name/bin/activate
```
3. install required python packages within the virtual environment
```bash
pip install -r path_to_this_dir/EnvironmentSetup/requirements_auto_process.txt
```
3. Locate your virtual environement and modify _config.sh_ by adding a line below
```
echo "export PYTHON_ENV='absolute_path_to_you_venv1_name/bin/activate'" > path_to_this_dir/config.sh
```

- **Usage** 
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
