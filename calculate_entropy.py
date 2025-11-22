#%% 

import os 
import numpy as np

import scipy

import pickle

from tqdm import tqdm


from scipy.special import lambertw
from scipy.stats import poisson 

#%%

from utils import return_none



#%%
############################################# 
## Calculate entropy for each cell barcode >>
#############################################
# @not-in-use 
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

    Mtotal = tn5_insert_array.shape[1] # total number of windows

    row_indices, col_indices, non0_insert =scipy.sparse.find(tn5_insert_array)
    Mnon0 = non0_insert.shape[0]


    non0_insert_adjusted = (non0_insert + 1) // 2 # adjust for paired-end sequencing
    cur_adj, freq_adj = np.unique(non0_insert_adjusted, return_counts=True)

    ZTP_prob = freq_adj / Mnon0
    mean_adj = np.matmul(cur_adj, ZTP_prob) 

    if mean_adj <= 1:
        return None, None, Mnon0, None, None, None

    # Lambert W function requires (-M1 * exp(-M1)) > -1/e --> mean_adj > 1 
    mle_lambda = lambertw( -mean_adj * np.exp(-mean_adj) , k=0).real + mean_adj
    P0 = np.exp(-mle_lambda)

    M0_openregion = Mnon0 / ( 1- P0) * P0
    p_closestate = 1 -  (M0_openregion + Mnon0) / Mtotal

    if p_closestate < 0:
        return None, None, Mnon0, mle_lambda, P0, p_closestate

    # Entropy of mix states
    p_openstate = 1 - p_closestate
    Entropy_states = -np.log2( p_closestate ) * p_closestate - p_openstate * np.log2( p_openstate )

    # Entropy only for open regions under Poisson model
    Poisson_prob = ZTP_prob * (1 - P0) # Only for non-zero windows in open region
    Entropy_openregion = -np.matmul(Poisson_prob, np.log2(Poisson_prob)) - P0 * np.log2(P0) 

    #Entropy_mixturedist = Entropy_states + p_openstate * Entropy_openregion
    Entropy_mixturedist = Entropy_states + p_openstate * Entropy_openregion / Mnon0 * (1 - P0) # Adjust for sample size effect


    return Entropy_mixturedist, Entropy_openregion, Mnon0, mle_lambda, P0, p_closestate



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
    barcode_entropy_df_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy_df.tsv')


    #%%
    ## Multiprocessing implementation
    if os.cpu_count() and os.cpu_count() >4:
        num_cpus = 4

        # Wrapper function for multiprocessing
        def mp_wrapper(bc):
            Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = mixdist_mle_entropy(insert_record[bc])
            return ( bc, Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state )

        from multiprocessing import Pool
        pool = Pool(processes=num_cpus)         # start num_cpus worker processes
        keys = list(insert_record.keys())
        res = pool.imap_unordered(mp_wrapper, keys, chunksize=4000)

        with open(barcode_entropy_df_file, 'w') as f:
            f.write(",Entropy,Entropy_open_region,Mp,mle_lambda,P0_open_region,P_closed_state\n")
            for r in tqdm(res, total=len(keys), desc="Calculating entropy for each cell barcode"):
                bc, Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = r
                # if np.isnan(Entropy_mixturedist):
                #     continue
                barcode_entropy[bc] = Entropy_mixturedist
                f.write(f"{bc},{Entropy_mixturedist},{Entropy_open_region},{Mnon0},{mle_lambda},{p0},{p_closed_state}\n")


    ## Single-threaded implementation 
    else:
        with open(barcode_entropy_df_file, 'w') as f:
            f.write(",Entropy,Entropy_open_region,Mp,mle_lambda,P0_open_region,P_closed_state\n")
            for bc, v in tqdm(insert_record.items(), desc="Calculating entropy for each cell barcode"):
                Entropy_mixturedist, Entropy_open_region, Mp, mle_lambda, p0, p_closed_state = mixdist_mle_entropy(v)
                # if np.isnan(Entropy_mixturedist):
                #     continue
                barcode_entropy[bc] = Entropy_mixturedist
                f.write(f"{bc},{Entropy_mixturedist},{Entropy_open_region},{Mp},{mle_lambda},{p0},{p_closed_state}\n")


    # save entropy results
    barcode_entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')
    with open(barcode_entropy_file, 'wb') as file:
        pickle.dump(barcode_entropy, file, protocol=pickle.HIGHEST_PROTOCOL)