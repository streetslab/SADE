#%% 

import os
import re 
import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix
from collections import defaultdict
import pickle
from tqdm import tqdm


## Global Variables >>
WindowSize = 3000  # size of the window in bp, default is 250 bp
## Global Variables <<

# Use the current directory for imports
from utils import return_none, get_chromsize_dict, assign_window, default_chrome_bin_sparse


# %%


# %%
def insert_frequency(chromosome:str, 
                     bl_region_file:str,
                     fragments_file:str, 
                     chromsize_dict:dict, 
                     output_file:str,
                     window_size:int = 3000):
    """Save insertion frequency for each cell barcode in the given chromosome.

    Args:
        fragments_file (str): tsv file 
        record_file (str): 
    """
    
    if chromosome +'_' not in fragments_file:
        raise ValueError(f"fragments_file: **{fragments_file}** does not contain chromosome: **{chromosome}**")
    
    bl = pd.read_csv(bl_region_file, sep='\t', header=None, usecols=[0, 1, 2], comment="#")
    bl.columns = ['chrom', 'start', 'end']

    windowdrop_indices = []
    for s, e in bl.loc[bl['chrom'] == chromosome, ['start', 'end']].itertuples(index=False):
        windowdrop_indices.append((s // window_size,  (e - 1) // window_size + 1))
    
    record = defaultdict(return_none)
    
    bc_list = {}
    with open(fragments_file, "r") as f:
        for l in tqdm(f, desc="Counting Tn5 insertion frequency"):
            if l.startswith("#"):
                continue
            chrom, left_insert, right_insert, cb, _ = l.strip().split("\t")  # cb -- cell barcode
            
            left_insert = int(left_insert)  # start-position is inclusive
            right_insert = int(right_insert) - 1  # end-position is exclusive

            if record[cb] is None:
                record[cb] = default_chrome_bin_sparse(chromsize_dict[chrom], window_size) # initialize the dictionary with zeros

            ind1 = assign_window(left_insert, chromsize_dict[chrom], window_size)
            ind2 = assign_window(right_insert, chromsize_dict[chrom], window_size)
            record[cb][0, ind1] += 1
            record[cb][0, ind2] += 1
            
            bc_list[cb] = True
        
    # Remove blacklisted windows from the record
    if windowdrop_indices:
        for cb in tqdm(bc_list, desc=f'Subtracting insertions in blacklisted regions from the record for chromosome {chromosome} ....'):
            for start, end in windowdrop_indices:
                record[cb][0, start:end] = 0
    else:
        print(f'All insertions are counted as no blacklisted regions provided for this chromosome {chromosome}.')
    
    # Save the insertion frequency record to a pickle file
    with open(output_file, 'wb') as f:
        pickle.dump(record, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    del record

# %%

if __name__ == "__main__":
    '''
    running example:
    python3 calculate_entropy.py --output_dir 'res' \
        --genome_chromsize 'ref/human_genome_chromsize.tsv' --chromosome 'chr22' 
    
    # input arguments:
    # output_dir: dir to save entropy calculation results 
    # frag_file: tsv file of fragments from CR output 
    # genome_chromsize: tsv file of human genome chromosome size
    
    '''

    import argparse
    
    parser = argparse.ArgumentParser(description="Calculate entropy for cell barcodes based on fragment file.\n  Entropy calculation is for one chromosome at a time.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save entropy calculation results.")
    parser.add_argument("--genome_chromsize_file", type=str, required=True, help="TSV file of human/species genome chromosome size.")
    parser.add_argument("--blacklist_file", type=str, required=True, help="bed file of blacklisted regions.")
    parser.add_argument("--chromosome", type=str, default="chr1", help="Chromosome to analyze, default is chr1.")
    parser.add_argument("--windowsize", type=int, default=WindowSize, help="Window size for fragmentation.")
    parser.add_argument("--shortcut_mode", type=str, required=False,  default="no", help="Go straight to Tn5 insertion counting and skip fragment splitting, default is no. no or yes")

    args = parser.parse_args()
    
    output_dir = args.output_dir
    genome_chromsize_file = args.genome_chromsize_file
    blacklist_file = args.blacklist_file
    windowsize = args.windowsize
    chromosomes_spec = args.chromosome

    
    chromsize_dict = get_chromsize_dict(genome_chromsize_file)
    
    shortcut_mode = args.shortcut_mode.lower() == "yes"
    if not shortcut_mode:
        save_subdir = output_dir
    else:
        save_subdir = os.path.join(output_dir, 'per_chromosome_entropy')
        if not os.path.exists(save_subdir):
            os.makedirs(save_subdir, exist_ok=True)
    chromosome_fragment_dir = os.path.join(output_dir, 'chromosome_fragments')
    
    
    # Do a file sanity check before counting even starts        
    if chromosomes_spec.lower() == "autosomes":
        chromosomes = []
        for chrom in chromsize_dict.keys():
            if re.match(r"^chr[0-9]+$", chrom) is not None:
                chromosomes.append(chrom)
    else:
        chromosomes = chromosomes_spec.split(' ')
        
    for chromosome in chromosomes:
        frag_file = os.path.join(chromosome_fragment_dir, f'{chromosome}_fragments.tsv')
        if not os.path.exists(frag_file):
            raise FileNotFoundError(f"Fragment file not found: {frag_file}")
    
    for chromosome in chromosomes:
        frag_file = os.path.join(chromosome_fragment_dir, f'{chromosome}_fragments.tsv')
            
        print(f"Counting Tn5 insertions for \n \
                    {chromosome} \n \
                using fragments from \n \
                    {frag_file} \n \
                and saving results to \n \
                    {save_subdir} .......... \n")
        
        print(f'chromosome size for {chromosome}: {chromsize_dict[chromosome]} \n ')
        insert_frequency_file = os.path.join(save_subdir, f'Pre_blacklist_{chromosome}_insert_frequency.pickle')
        
        # Step 2. Count Tn5 insertion frequency for each cell barcode in the given chromosome
        insert_frequency(chromosome=chromosome,
                         bl_region_file=blacklist_file,
                        fragments_file=frag_file,
                        chromsize_dict=chromsize_dict,
                        output_file=insert_frequency_file,
                        window_size=windowsize)
        


