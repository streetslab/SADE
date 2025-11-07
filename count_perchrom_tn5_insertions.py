#%% 

import os 
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
from utils import return_none


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
def insert_frequency(chromosome:str, 
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
    
    
    record = defaultdict(return_none)
    
    
    with open(fragments_file, "r") as f:
        #print(f"counting tn5 insertion frequency for {chromosome} in {os.path.basename(fragments_file)} ..........")
        for l in tqdm(f, desc="Counting tn5 insertion frequency"):
            if l.startswith("#"):
                continue
            chrom, left_insert, right_insert, cb, _ = l.strip().split("\t")  # cb -- cell barcode
            

            left_insert = int(left_insert)
            right_insert = int(right_insert)

            if record[cb] is None:
                total_seq = int(chromsize_dict[chrom] // window_size + 1)
                record[cb] = lil_matrix((1, total_seq), dtype=np.int32) # initialize the dictionary with zeros

            ind1 = assign_window(left_insert, chromsize_dict[chrom], window_size)
            ind2 = assign_window(right_insert, chromsize_dict[chrom], window_size)
            record[cb][0, ind1] += 1
            record[cb][0, ind2] += 1
    
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
    parser.add_argument("--genome_chromosize_file", type=str, required=True, help="TSV file of human/species genome chromosome size.")
    parser.add_argument("--chromosome", type=str, default="chr1", help="Chromosome to analyze, default is chr1.")
    parser.add_argument("--windowsize", type=int, default=WindowSize, help="Window size for fragmentation.")

    args = parser.parse_args()
    
    output_dir = args.output_dir
    genome_chromsize_file = args.genome_chromosize_file
    chromosome = args.chromosome
    windowsize = args.windowsize
    
    chromosome_fragment_dir = os.path.join(output_dir, 'chromosome_fragments')
    frag_file = os.path.join(chromosome_fragment_dir, f'{chromosome}_fragments.tsv')

    if not os.path.exists(frag_file):
        raise FileNotFoundError(f"Fragment file not found: {frag_file}")

    print(f"Calculating entropy for \n \
                {chromosome} \n \
            using fragments from \n \
                {frag_file} \n \
            and saving results to \n \
                {output_dir} .......... \n")
    

    # load human genome chromosome size
    chromsize_dict = {}
    with open(genome_chromsize_file, "r") as f:
        for line in f:
            chrom, _, size = line.strip().split("\t")
            chromsize_dict[chrom] = int(size)

    print(f'chromosome size for {chromosome}: {chromsize_dict[chromosome]} \n ')
    insert_frequency_file = os.path.join(output_dir, f'{chromosome}_insert_frequency.pickle')
    

    # Step 2. Count Tn5 insertion frequency for each cell barcode in the given chromosome
    insert_frequency(chromosome=chromosome,
                     fragments_file=frag_file,
                     chromsize_dict=chromsize_dict,
                     output_file=insert_frequency_file,
                     window_size=windowsize)


