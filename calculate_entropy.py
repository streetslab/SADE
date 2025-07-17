#%% 

import os 
import numpy as np
import pandas as pd
import psutil 
from scipy.sparse import lil_matrix
from collections import defaultdict
import gc
import copy
import pickle
import matplotlib.pyplot as plt
from tqdm import tqdm

## Global Variables >>
WindowSize = 250  # size of the window in bp, default is 250 bp
#precision = 1e-10  # to avoid log2(0) in entropy calculation
precision = 0
## Global Variables <<


# %%


def assign_window(insertion_pos, chromsize, window_size):
    """
    Assign tn5_insertion position to a window
    Args:
        insertion_pos (int): tn5_insertion position
        chromsize (int): chromosome size
        window_size (int): window size
    Returns: index of the window that needs to be incremented with 1 
    """
    if insertion_pos < 0 or insertion_pos >= chromsize:
        raise ValueError(f"Invalid insertion position: {insertion_pos} for chromosome size: {chromsize}")
    
    return insertion_pos // window_size 

# %%
def insert_frequency(fragments_file:str, 
                     record = None, 
                     chromsize_dict:dict = None,
                     chromosome:str = "chr1",
                     window_size:int = 250):
    """return insertion frequency for each cell barcode in the given chromosome.
    (1) If record is None, it will return a new record with cell barcodes as keys and insertion frequency as values.
    The record is a dictionary with cell barcodes as keys and sparse matrix as values.
    
    (2) If record is not None, it will update the record in place.

    Args:
        fragments_file (str): tsv file 
        record (dict): None of defaultdict
    """
    
    return_record = False
    if record is None:
        return_record = True
        record = defaultdict(lambda: None)
    with open(fragments_file, "r") as f:
        print(f"counting tn5 insertion frequency for {chromosome} in {fragments_file}")
        for l in tqdm(f, desc="Counting tn5 insertion frequency"):
            if l.startswith("#"):
                continue
            chrom, left_insert, right_insert, cb, _ = l.strip().split("\t")  # cb -- cell barcode
            
            if chrom == chromosome:
                left_insert = int(left_insert)
                right_insert = int(right_insert)

                if record[cb] is None:
                    total_seq = int(chromsize_dict[chrom] // window_size + 1)
                    record[cb] = lil_matrix((1, total_seq), dtype=np.int32) # initialize the dictionary with zeros

                ind1 = assign_window(left_insert, chromsize_dict[chrom], window_size)
                ind2 = assign_window(right_insert, chromsize_dict[chrom], window_size)
                record[cb][0, ind1] += 1
                record[cb][0, ind2] += 1
            
    if return_record:
        return record
    else:
        return None


# %%

if __name__ == "__main__":
    '''
    running example:
    python3 calculate_entropy.py --res_dir 'res' \
        --frag_file '/home/syyang/adipose_ln/atac/_cr_atac_results/eda_AB_ATAC_FL_TAM_SQ/fragments.tsv'  \
        --genome_chromsize 'ref/human_genome_chromsize.tsv' --chromosome 'chr22'
    
    # input arguments:
    # res_dir: dir to save entropy calculation results 
    # frag_file: tsv file of fragments from CR output 
    # genome_chromsize: tsv file of human genome chromosome size
    
    '''

    import argparse
    
    parser = argparse.ArgumentParser(description="Calculate entropy for cell barcodes based on fragment file.\n  Entropy calculation is for one chromosome at a time.")
    parser.add_argument("--res_dir", type=str, required=True, help="Directory to save entropy calculation results.")
    parser.add_argument("--frag_file", type=str, required=True, help="TSV file of fragments from CR output.")
    parser.add_argument("--genome_chromsize", type=str, required=True, help="TSV file of human/species genome chromosome size.")
    parser.add_argument("--chromosome", type=str, default="chr1", help="Chromosome to analyze, default is chr1.")

    args = parser.parse_args()
    
    res_dir = args.res_dir
    frag_file = args.frag_file
    genome_chromsize_file = args.genome_chromsize
    chromosome = args.chromosome
    
    
    print(f"Calculating entropy for \n \
            {chromosome} \n \
            using fragments from \n \
            {frag_file} \n \
            and saving results to \n \
                {res_dir}")

    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    
    import subprocess
    # Step 1. Filter fragments based on chromosome in shell script 
    chrome_frag_file = os.path.join(res_dir, f"{chromosome}_fragments.tsv")
    subprocess.run(f"grep ^{chromosome} {frag_file} > {chrome_frag_file}", shell=True)
    
    # TODO: check if the chrome_frag_file is empty, if so, raise an error. --stop here. 
    
    
    # load human genome chromosome size
    chromsize_dict = {}
    with open(genome_chromsize_file, "r") as f:
        for line in f:
            chrom, _, size = line.strip().split("\t")
            chromsize_dict[chrom] = int(size)

    print(f'chromosome size for {chromosome}: {chromsize_dict[chromosome]}')
    
    

    # Step 2. Count Tn5 insertion frequency for each cell barcode in the given chromosome
    _insert_record = defaultdict(lambda: None)
    _  = insert_frequency(chrome_frag_file, _insert_record, chromsize_dict, chromosome)  # inplace update of the record -- save memory. IDK why returning record is not memory efficient.
    
    # convert defaultdict to dict  and save result
    insert_record = {k: v for k, v in _insert_record.items() } 

    insert_frequency_file = os.path.join(res_dir, f'{chromosome}_insert_frequency.pickle')
    with open(insert_frequency_file, 'wb') as file:
        pickle.dump(insert_record, file, protocol=pickle.HIGHEST_PROTOCOL)
        

    # %%
    # Step 3. Calculate entropy for each cell barcode for the given chromosome 
    barcode_entropy = {}
    for k, v in tqdm(insert_record.items(), desc="Calculating entropy for each cell barcode"):
        occur, freq = np.unique(v.toarray(), return_counts=True)
        freq = freq / np.sum(freq)
        entropy = -np.sum(freq * np.log2(freq + precision))
        barcode_entropy[k] = entropy

    # save entropy results
    barcode_entropy_file = os.path.join(res_dir, f'{chromosome}_barcode_entropy.pickle')
    with open(barcode_entropy_file, 'wb') as file:
        pickle.dump(barcode_entropy, file, protocol=pickle.HIGHEST_PROTOCOL)
        
    # Step 4 (optional -- visualization). Make Knee plot of entropies 
    figure_subdir = os.path.join(res_dir, 'figures')
    if not os.path.exists(figure_subdir):
        os.makedirs(figure_subdir)
    import matplotlib.pyplot as plt
    
    entropies = np.array(list(barcode_entropy.values()))
    entropies.sort() # sort in place and ascending order 
    fig, ax = plt.subplots()
    ax.plot(entropies[::-1], marker='o', linestyle='-', markersize=2)
    ax.set_xlabel('Cell Barcode Index (sorted by entropy)')
    ax.set_ylabel('Entropy')
    fig.suptitle(f'Entropy distribution for {chromosome}')
    fig.savefig(os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot.png'))
    plt.close(fig)  