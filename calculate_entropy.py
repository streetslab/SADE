#%% 

import os 
import numpy as np
import pandas as pd
import psutil 
import scipy
from scipy.sparse import lil_matrix
from collections import defaultdict
import gc
import copy
import pickle
import matplotlib.pyplot as plt
from tqdm import tqdm


from scipy.special import lambertw
from scipy.stats import poisson 

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
    parser.add_argument("--windowsize", type=int, default=WindowSize, help="Window size for fragmentation, default is 0.001.")

    args = parser.parse_args()
    
    res_dir = args.res_dir
    frag_file = args.frag_file
    genome_chromsize_file = args.genome_chromsize
    chromosome = args.chromosome
    windowsize = args.windowsize
    
    
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
    _  = insert_frequency(chrome_frag_file, _insert_record, chromsize_dict, chromosome, windowsize)  # inplace update of the record -- save memory. IDK why returning record is not memory efficient.
    
    # convert defaultdict to dict  and save result
    insert_record = {k: v for k, v in _insert_record.items() } 

    insert_frequency_file = os.path.join(res_dir, f'{chromosome}_insert_frequency.pickle')
    with open(insert_frequency_file, 'wb') as file:
        pickle.dump(insert_record, file, protocol=pickle.HIGHEST_PROTOCOL)


    ############################################# 
    ## Calculate entropy for each cell barcode >>
    #############################################
    #%% 
    ## MLE estimation of Poisson distribution parameter with Zero-Truncation Poisson observations
    def mle_estimate(tn5_insert_array:scipy.sparse, print_debug:bool=False, return_freq:bool=False):
        flag = None
        Mtotal = tn5_insert_array.shape[1] # total number of windows
        
        row_indices, col_indices =tn5_insert_array.nonzero()
        non0_insert = tn5_insert_array[row_indices, col_indices].toarray().flatten()
        Mp = non0_insert.shape[0]
        
        cur, freq = np.unique(non0_insert, return_counts=True)
        
        non0_insert_adjusted = (non0_insert + 1) // 2 # adjust for paired-end sequencing
        cur_adjusted, freq_adjusted = np.unique(non0_insert_adjusted, return_counts=True)
        
        p_non0 = freq_adjusted / Mp
        mean_adjusted = np.sum( cur_adjusted * p_non0 )

        # Lambert W function requires (-M1 * exp(-M1)) >= -1/e
        if mean_adjusted <= 1:
            mle_lambda = None
            P0 = None
            if return_freq:
                return Mtotal, Mp, mean_adjusted, mle_lambda, P0, cur_adjusted, freq_adjusted, cur, freq, flag
            return Mtotal, Mp, mean_adjusted, mle_lambda, P0, cur_adjusted, freq_adjusted, cur, freq, flag

        mle_lambda = lambertw( -mean_adjusted * np.exp(-mean_adjusted) , k=0).real + mean_adjusted
        P0 = np.exp(-mle_lambda)
        
        if mle_lambda <= 0:
            P0 = None
        
        if return_freq:
            return Mtotal, Mp, mean_adjusted, mle_lambda, P0, cur_adjusted, freq_adjusted, cur, freq, flag
        return Mtotal, Mp, mean_adjusted, mle_lambda, P0, cur_adjusted, freq_adjusted, cur, freq, flag


    def kl_divergence(mle_lambda, cur, freq):
        """ Compute KL divergence D_KL(P||Q) for discrete distributions
            p and q are arrays of the same length representing probability distributions
        """
        # q is from Poisson distribution with parameter mle_lambda
        
        p = freq / np.sum(freq)
        log_q = np.array([poisson.logpmf(k, mle_lambda) for k in cur]) 
        
        divergence = np.sum( p * (np.log(p) - log_q) )
        
        return divergence
    

    # Calculate Entropy using MLE estimate 
    def mixdist_mle_entropy(tn5_insert_array:scipy.sparse):
        Mtotal, Mp, mean_non0, mle_lambda, p0, cur_adjusted, freq_adjusted,  cur, freq, flag = mle_estimate(tn5_insert_array, return_freq=True)

        if p0 is None:
            Entropy_mixturedist = None
            Entropy_open_region = None
            Entropy_open_region_noadj = None
            divergence = None
            p_closed_state = None
        else:
            Nwindows_0inser_open_region = np.floor(Mp / ( 1- p0) * p0 )
            Nwindows_0inser_open_region = np.min((Nwindows_0inser_open_region, Mtotal - Mp))
            Nwindows_closed_region = Mtotal - Mp - Nwindows_0inser_open_region

            #total_freq = np.concatenate(([Nwindows_closed_region, Nwindows_0inser_open_region], freq_adjusted))

            # calculate entropy of mix states
            p_closed_state = Nwindows_closed_region / Mtotal
            p_states = np.array([p_closed_state, 1 - p_closed_state])
            Entropy_state = -np.sum(p_states * np.log2(p_states))
            
            # calculate entropy as a mixture distribution 
            if Nwindows_0inser_open_region > 0:
                open_region_p = np.concatenate( ([Nwindows_0inser_open_region ], freq_adjusted) )
            else:
                open_region_p = freq_adjusted
            open_region_p = open_region_p / (Mp + Nwindows_0inser_open_region)
            Entropy_open_region = -np.sum(open_region_p * np.log2(open_region_p))

            # # Calculate the open region entropy using unadjusted freq (for comparison)
            p = freq/Mp
            Entropy_open_region_noadj = -np.sum(p * np.log2(p))
            
            Entropy_mixturedist = Entropy_state + (1 - p_closed_state) * Entropy_open_region
            # flag =  Entropy_state + (1 - p_closed_state) * Entropy_open_region_noadj # to test the difference
            #Entropy_mixturedist = Entropy_state + (1 - p_closed_state) * entropy_non0_part # without adjusted
            divergence = kl_divergence(mle_lambda, cur_adjusted, freq_adjusted)
        return Entropy_mixturedist, Entropy_open_region, Entropy_open_region_noadj, mle_lambda, p0, p_closed_state, divergence


    ############################################# 
    ## Calculate entropy for each cell barcode <<
    #############################################

    # %%
    # Step 3. Calculate entropy for each cell barcode for the given chromosome 
    barcode_entropy = {}
    entropy_df = pd.DataFrame( columns=[
        'entropy_open_region', 
        'Entropy_open_region_noadj', 
        'mle_lambda', 
        'P0_open_region', 
        'P_closed_state',
        'kl_divergence'
    ])
    
    for bc, v in tqdm(insert_record.items(), desc="Calculating entropy for each cell barcode"):
        Entropy_mixturedist, Entropy_open_region, Entropy_open_region_noadj, mle_lambda, p0, p_closed_state, divergence = mixdist_mle_entropy(v)
        barcode_entropy[bc] = Entropy_mixturedist
        entropy_df.loc[bc, 'Entropy'] = Entropy_mixturedist
        entropy_df.loc[bc, 'Entropy_open_region'] = Entropy_open_region
        entropy_df.loc[bc, 'Entropy_open_region_noadj'] = Entropy_open_region_noadj
        entropy_df.loc[bc, 'mle_lambda'] = mle_lambda
        entropy_df.loc[bc, 'P0_open_region'] = p0
        entropy_df.loc[bc, 'P_closed_state'] = p_closed_state
        entropy_df.loc[bc, 'kl_divergence'] = divergence


    # save entropy results
    barcode_entropy_file = os.path.join(res_dir, f'{chromosome}_barcode_entropy.pickle')
    with open(barcode_entropy_file, 'wb') as file:
        pickle.dump(barcode_entropy, file, protocol=pickle.HIGHEST_PROTOCOL)
        
    barcode_entropy_df_file = os.path.join(res_dir, f'{chromosome}_barcode_entropy_df.tsv')
    entropy_df.to_csv(barcode_entropy_df_file, sep='\t')
        
    # Step 4 (optional -- visualization). 
    ## Make Knee plot of entropies 
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
    
    
    ## Make Knee plot of entropies on log10 scale
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.log10(entropies[::-1] + precision), '-', color='blue', alpha=0.5, label='Entropy')
    ax.set_ylabel('Entropy (log10-scaled)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_logscale.png')
    fig.savefig(fig_file, bbox_inches='tight')
    

    #histogram of log10(entropy) values -- for (potentially) Gaussian mixture model fitting
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(np.log10(entropies + precision), bins=100, color='blue',  label='log10(Entropy)')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('log10(Entropy)')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_histogram.png')
    fig.savefig(fig_file, bbox_inches='tight')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(np.log10(entropies + precision), bins=100, color='blue',  log=True, label='log10(Entropy)')
    ax.set_ylabel('log10(Frequency)')
    ax.set_xlabel('log10(Entropy)')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_histogram_logFreq.png')
    fig.savefig(fig_file, bbox_inches='tight')