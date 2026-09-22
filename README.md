# SADE
**S**hannon-entropy for sc**A**TAC-seq **De**noising.

**SADE** quantifies the complexity of genome-wide chromatin accessibility for each droplet(-barcode) with Shannon entropy, and auto-threshold entropy values to identify droplets with quality nuclei, while filtering out empty droplets, droplets with damaged nuclei, as well as droplets with cellular-debris.  


## System Requirements
[ripgrep](https://github.com/BurntSushi/ripgrep)   
[python3](https://www.python.org/downloads/)  


## Runtime
SADE runs in a few minutes on standard macOS and Linux systems with 16 GB RAM, depending on dataset size and selected parameters.     
As reference: with default genome window-size at 3000, an experiment that targeted for 10k human cells took 20-min to finish running SADE, an experiment that targeted for 2k human cells took 5-min to finish running SADE.    
Specie's genome size is positively correlated with run-time, specified genome binning window-size is negatively correlated with run-time.     
If you don't see progress messages popping up soon after starting running SADE, it is not properly set up. 

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
          ```-c <chromosome>``` (Optional): Chromosome used to calculate entropy, space-separated list, or a single chromosome, e.g. 'chr1', 'chr1 chr2'. Default: chr1.  
          ```-w <window_size>``` (Optional): Window size (in basepair) on genome to count Tn5 insertion frequencies. Default: 3000.    
                **Note:** Window-size 3000 is recommended for human/mouse genomes, although sizes between 500-10000 were tested as robust for these two species. We recommend roughly scale window size to genome size for other species.    
          ```-k <candidate_inflection_points>```  (Optional): Number of candidate local inflection points being considered to find the global optimal inflection point on entropy value fitted cubic curve. Default: 2.  
          ```-s <genome_saturation_cutoff>``` (Optional): Upper bound of the genome portion that can possibly have experimental (i.e. Tn5 insertion) signals for any cells. Default: 0.5.  

  * **Output**  
      - Main:  
          ```fragments.tsv```:   A copy of the input fragment file.  
          ```filtered_fragments.tsv```:  Fragments corresponding to quality nuclei based on the entropy criterion. Use this file for downstream analysis.  
          ```entropy_filtered_bc_df.tsv```:  A .tsv (tab-separated) file containing filtered barcodes as row indices and entropy calculation metrics as columns.  
          ```figures```: A subfolder containing generated plots.  
      - Additional:   
          ```Entropy_filtered_bc_CBZ.txt```: A text file with each row being SADE identified high quality nuclei (droplet-)barcode with added CBZ tag (e.g. "CB:Z:ATTTGCAAGTATTGTG-1"). All the scATACseq datasets used are from 10x, so 10x specified tag pre-pending the filtered barcodes are saved in a file. This file can then be used to filter bam/sam file if needed. Although the filtered barcodes are in the ```entropy_filtered_bc_df.tsv``` and can be processed in alternative ways such as adding different tags to be compatible for other protocols.
        <ins>**TODO for myself to make it more general/compatible to other protocols on this file?**</ins>
    
  * **Note**   
            Currently 'hg38' and 'mm10' genome reference info are pre-processed and stored in the ```ref``` subfolder in this repo.  
            If cells come from a different species or a different genome version needs to be used, users can refer to ```ref/get_chromosome_size.sh``` to get genome size information and store it with the same file naming convention (i.e. ```{new_species_name}_genome_chromsize.tsv```) into the ```ref/``` subfolder in your local path.  
            You can then use ```sade.sh``` command demonstrated above to specify ```-g {new_species_name}``` 
    
## Filter bam file (Optional)
Filter .bam file that corresponds to the filtered fragment file at (droplet-)barcode level (if needed)  
[samtools](https://www.htslib.org) is required and callable in the system for this step.  
* **Note**   
Users will need the output file ```Entropy_filtered_bc_CBZ.txt``` from running ```sade.sh``` to filter bam file.   
* **Note**   
```filter_bam_with_barcodes.sh``` file here is for filtering scATACseq bam from 10x protocols given its unique barcode tag, 
and code in this shell script was adapted from [10x](https://kb.10xgenomics.com/s/article/360022448251-How-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes).    
Filtering on bam files generated from other experimental protocols should be able to use this script as reference.

```bash
filtered_bc_file="Entropy_filtered_bc_CBZ.txt" # See SADE's output above.
bash path_to_this_dir/filter_bam_with_barcodes.sh \
            -s ${unfiltered_bam_file} \
            -o ${you_specified_output_dir} \
            -f ${filtered_bc_file}
```
