# SADE
**S**h**a**nnon-entropy **de**tection of quality nuclei in scATAC-seq.

**SADE** quantifies the complexity of genome-wide chromatin accessibility for each droplet(-barcode) with entropy, and auto-threshold entropy values to identify droplets with quality nuclei, while filtering out empty droplets, droplets with damaged nuclei, as well as droplets with cellular-debris.  


## System Requirements
[ripgrep](https://github.com/BurntSushi/ripgrep)


## Environment Setup

To run ```sade.sh```, a Python 3 virtual environment is required to calculate per-(droplet-)barcode entropy and perform auto-thresholding.
1. Create a python virtual environement with required libraries 
```bash 
python -m venv you_venv1_name
```
2. Activate your python virtual environment and install required python packages within the virtual environment
```bash
source absolute_path_to_you_venv1_name/bin/activate
pip install -r path_to_this_dir/EnvironmentSetup/pip_requirements.txt
```
3. Locate your virtual environement and modify _config.sh_ by adding a line below
```
echo "export PYTHON_ENV='absolute_path_to_you_venv1_name/bin/activate'" > path_to_this_dir/config.sh
```

## Usage 
```
sample_dir='path_to_your_sample_dir' # Need to contain fragment file (fragments.tsv or fragments.tsv.gz).
ws='3000'
bash path_to_this_dir/sade.sh -o your_desired_output_directory \
      -f "${sample_dir}/fragments.tsv.gz" \
      -w "${ws}" \
      -c 'chr1' \
      -g 'hg38' # hg38 for the corresponding human sample
  ```   
  
### Command **Explanation**:
  * **Input**  
          ```-o <output_dir>``` (Required): Path to the output directory.  
          ```-f <fragments_file>``` (Required): Path to the fragments file.  
          ```-g <genome>```  (Required): Genome used for read mapping (e.g., hg38, mm10).  
          ```-c <chromosome>``` (Optional): Chromosome used to calculate entropy. Default: chr1.  
          ```-w <window_size>``` (Optional): Window size (in basepair) on genome to count Tn5 insertion frequencies. Default: 3000.  
          ```-k <candidate_inflection_points>```  (Optional): Number of candidate local inflection points being considered to find the global optimal inflection point on entropy value fitted cubic curve. Default: 2.  
          ```-s <genome_saturation_cutoff>``` (Optional): Upper bound of the genome portion that can possibly have experimental (i.e. Tn5 insertion) signals for any cells. Default: 0.5.  

  * **Output**   
          ```fragments.tsv```:   A copy of the input fragment file.  
          ```filtered_fragments.tsv```:  Fragments corresponding to quality nuclei based on the entropy criterion. Use this file for downstream analysis.  
          ```entropy_filtered_bc_df.tsv```:  A .tsv (tab-separated) file containing filtered barcodes as row indices and entropy calculation metrics as columns.  
          ```figures```: A subfolder containing generated plots.  
    
  * **Note**   
            Currently 'hg38' and 'mm10' genome reference info are pre-processed and stored in the ```ref``` subfolder in this repo.  
            If cells come from a different species or a different genome version needs to be used, users can refer to ```ref/get_chromosome_size.sh``` to get genome size information and store it in the same file naming convention (i.e. ```{new_species_name}_genome_chromsize.tsv```) approach to the ```ref/``` subfolder in your local path.  
            You can then use **SADE** command demonstrated above to specify ```-g {new_species_name}``` 
    
## Filter bam file (Optional)
Filter .bam file that corresponds to the filtered fragment file at (droplet-)barcode level (if needed)  
[samtools](https://www.htslib.org) is required for this step.

**TODO**
**ADD sh files back for this**
