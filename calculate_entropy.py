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

precision = 1e-7 # avoid nan in lambert W function calculation

# %%

# For loading pickle files of Tn5 insertion frequency that has defaultdict structure
from utils import return_none



#%%
############################################# 
## Calculate entropy for each cell barcode >>
#############################################

## MLE estimation of Poisson distribution parameter with Zero-Truncation Poisson observations
def mle_estimate(tn5_insert_array:scipy.sparse, print_debug:bool=False, return_freq:bool=False):
    
    Mtotal = tn5_insert_array.shape[1] # total number of windows
    
    row_indices, col_indices, non0_insert =scipy.sparse.find(tn5_insert_array)
    Mp = non0_insert.shape[0]
    
    
    non0_insert_adjusted = (non0_insert + 1) // 2 # adjust for paired-end sequencing
    cur_adjusted, freq_adjusted = np.unique(non0_insert_adjusted, return_counts=True)
    
    p_non0 = freq_adjusted / Mp
    mean_adjusted = np.sum( cur_adjusted * p_non0 ) + precision # make mean > 1 for lambertw function

    # Lambert W function requires (-M1 * exp(-M1)) > -1/e --> mean_adjusted > 1 
    mle_lambda = lambertw( -mean_adjusted * np.exp(-mean_adjusted) , k=0).real + mean_adjusted
    P0 = np.exp(-mle_lambda)
    
    if return_freq:
        return Mtotal, Mp, mle_lambda, P0, cur_adjusted, freq_adjusted
    return Mtotal, Mp, mle_lambda, P0, cur_adjusted, freq_adjusted


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
    Mtotal, Mp, mle_lambda, p0, cur_adjusted, freq_adjusted = mle_estimate(tn5_insert_array, return_freq=True)

    if p0 is None:
        Entropy_mixturedist = Entropy_open_region = divergence = p_closed_state = None
        return Entropy_mixturedist, Entropy_open_region, Mp, mle_lambda, p0, p_closed_state
    else:
        Nwindows_0inser_open_region = Mp / ( 1- p0) * p0 
        Nwindows_closed_region = Mtotal - Mp - Nwindows_0inser_open_region  # what if negative?

        if Nwindows_closed_region < 0:
            Entropy_mixturedist = Entropy_open_region = divergence = p_closed_state = None
            return Entropy_mixturedist, Entropy_open_region, Mp, mle_lambda, p0, p_closed_state
        
        
        # calculate entropy of mix states
        p_closed_state = Nwindows_closed_region / Mtotal
        Entropy_state = - p_closed_state * np.log2(p_closed_state) - (1 - p_closed_state) * np.log2(1 - p_closed_state)
        
        # calculate entropy as a mixture distribution 
        open_region_p = freq_adjusted / Mp * (1 - p0)
        Entropy_open_region = -np.matmul(open_region_p, np.log2(open_region_p)) - p0 * np.log2(p0)
        
        
        Entropy_mixturedist = Entropy_state + (1 - p_closed_state) * Entropy_open_region


        #divergence = kl_divergence(mle_lambda, cur_adjusted, freq_adjusted)
    return Entropy_mixturedist, Entropy_open_region, Mp, mle_lambda, p0, p_closed_state#, divergence



############################################# 
## Calculate entropy for each cell barcode <<
#############################################

# %%

if __name__ == "__main__":
    '''
    running example:
    python3 calculate_entropy.py --output_dir 'res' 
    
    # input arguments:
    # output_dir: dir to save entropy calculation results 
    # chromosome: chromosome to analyze, default is chr1    
    '''

    import argparse
    
    parser = argparse.ArgumentParser(description="Calculate entropy for cell barcodes based on fragment file.\n  Entropy calculation is for one chromosome at a time.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save entropy calculation results.")
    parser.add_argument("--chromosome", type=str, default="chr1", help="Chromosome to analyze, default is chr1.")

    args = parser.parse_args()
    
    output_dir = args.output_dir
    chromosome = args.chromosome
    
    

    # Step 1. Load Tn5 insertion frequency record for each cell barcode
    tn5_insert_record_file = os.path.join(output_dir, f'{chromosome}_insert_frequency.pickle')
    
    with open(tn5_insert_record_file, 'rb') as file:
        insert_record = pickle.load(file)

    # %%
    # Step 2. Calculate entropy for each cell barcode for the given chromosome 
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
        Entropy_mixturedist, Entropy_open_region, Mp, mle_lambda, p0, p_closed_state = mixdist_mle_entropy(v)
        barcode_entropy[bc] = Entropy_mixturedist
        entropy_df.loc[bc, 'Entropy'] = Entropy_mixturedist
        entropy_df.loc[bc, 'Entropy_open_region'] = Entropy_open_region
        entropy_df.loc[bc, 'Mp'] = Mp
        entropy_df.loc[bc, 'mle_lambda'] = mle_lambda
        entropy_df.loc[bc, 'P0_open_region'] = p0
        entropy_df.loc[bc, 'P_closed_state'] = p_closed_state


    # save entropy results
    barcode_entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')
    with open(barcode_entropy_file, 'wb') as file:
        pickle.dump(barcode_entropy, file, protocol=pickle.HIGHEST_PROTOCOL)
        
    barcode_entropy_df_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy_df.tsv')
    entropy_df.to_csv(barcode_entropy_df_file, sep='\t')
        
