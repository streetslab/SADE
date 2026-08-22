#%% 

import os 
import numpy as np

import scipy

import pickle

from tqdm import tqdm


from scipy.special import lambertw
from scipy.stats import poisson 
from scipy.sparse import hstack
from collections import defaultdict

import argparse

#%%

from utils import return_none



#%%
############################################# 
## Calculate entropy for each cell barcode >>
#############################################
def mixdist_mle_entropy(tn5_insert_array:scipy.sparse):
    # Use insertions >=2 to calculate poisson-rate for open regions. 
    Mtotal = tn5_insert_array.shape[1] # total number of windows

    row_indices, col_indices, non0_insert =scipy.sparse.find(tn5_insert_array)
    Mnon0 = non0_insert.shape[0]


    non0_insert_adjusted = (non0_insert + 1) // 2 # adjust for paired-end sequencing
    cur_adj, freq_adj = np.unique(non0_insert_adjusted, return_counts=True)
    
    # Use only insertions >=2 to calculate poisson-rate for open regions.
    if cur_adj[0] < 2:
        M1 = freq_adj[0]
        cur_adj_2 = cur_adj[1:]
        freq_adj_2 = freq_adj[1:]
    else:
        M1 = 0
        cur_adj_2 = cur_adj
        freq_adj_2 = freq_adj

    R_l2 = (cur_adj_2 * freq_adj_2).sum() # total number of insertions in open regions (across all windows with >=2 insertions)
    M_l2 = freq_adj_2.sum() # number of windows with >=2 insertions 

    def f_lambda(_lmbda, R_l2, M_l2):
        return (R_l2 - M_l2 * _lmbda) * np.exp(_lmbda) + ( 1+ M_l2) * _lmbda **2  - (2 - M_l2 - R_l2) * _lmbda  - R_l2
    
    #Newton-Raphson method can't gaurantee the non-negative solutions hence using Ridder's method to find the positive root
    mle_lambda = scipy.optimize.ridder(f_lambda, 1e-2, 100, (R_l2, M_l2))
    
    P_l2_openregion = 1 - (1 - mle_lambda) * np.exp(-mle_lambda) # P(X>=2) for Poisson distribution with rate mle_lambda
    M_est_openregion = M_l2 / P_l2_openregion # estimated number of total windows in open regions
    M1_est_openregion = M_est_openregion * mle_lambda *np.exp(-mle_lambda) # estimated number of windows with exactly 1 insertion in open regions
    M0_est_openregion = M1_est_openregion / mle_lambda # estimated number of windows with exactly 0 insertion in open regions

    # TODO: condition control 
    # If M1_est_openregion is >> M1 --> the model is not a good fit and is likely over-tagmentated dna
    if M1_est_openregion > M1:
        return None, None, None, None, mle_lambda, None, P_l2_openregion, Mnon0
    # If M0_est_openregion ???
    if M0_est_openregion > Mtotal - Mnon0:
        return None, None, None, None, mle_lambda, None, P_l2_openregion, Mnon0

    M0_est_closeregion = (Mtotal - Mnon0) - M0_est_openregion  
    M1_est_closeregion = M1 - M1_est_openregion
    lambda_closeregion = M1_est_closeregion / M0_est_closeregion # Poisson rate for closed regions (rate for background noise)

    if lambda_closeregion > mle_lambda:
        return None, None, None, None, mle_lambda, lambda_closeregion, P_l2_openregion, Mnon0
    # Use Pi_closestate for P(states=closed)
    Pi_closestate = 1 - M_est_openregion / Mtotal
    Pi_openstate = 1 - Pi_closestate
    Entropy_states = -np.log2( Pi_closestate ) * Pi_closestate - Pi_openstate * np.log2( Pi_openstate )

    # Entropy only for open regions under Poisson model
    # est_cur_adj = np.append([0, 1], cur_adj_2)
    est_freq_adj = np.append([M0_est_openregion, M1_est_openregion], freq_adj_2)
    Poisson_prob =  est_freq_adj / M_est_openregion  
    Entropy_openregion = -np.matmul(Poisson_prob, np.log2(Poisson_prob)) 
    
    # Entropy for closed region under Background Poisson model (small rate) 
    Poisson_prob_closeregion = np.array([M0_est_closeregion, M1_est_closeregion]) / (M0_est_closeregion + M1_est_closeregion)
    Entropy_closeregion =   -np.matmul(Poisson_prob_closeregion, np.log2(Poisson_prob_closeregion))

    Entropy_mixturedist = Entropy_states + Pi_closestate * Entropy_closeregion + Pi_openstate * Entropy_openregion / M_est_openregion # Adjust for sample size effect (indirectly read-depth)


    return Entropy_mixturedist, Entropy_openregion, Entropy_closeregion, Pi_closestate, mle_lambda, lambda_closeregion, P_l2_openregion, Mnon0

# Calculate Entropy using MLE estimate 
def mixdist_mle_entropy_ref(tn5_insert_array:scipy.sparse):

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
    Entropy_mixturedist = Entropy_states + p_openstate * Entropy_openregion / Mnon0 * (1 - P0) # Adjust for sample size effect (indirectly read-depth)


    return Entropy_mixturedist, Entropy_openregion, Mnon0, mle_lambda, P0, p_closestate


############################################# 
## Calculate entropy for each cell barcode <<
#############################################


##  >>
##  >> Move all the code before multiprocessing outside of __main__ 
parser = argparse.ArgumentParser(description="Calculate entropy for cell barcodes based on fragment file.\n  Entropy calculation is for one chromosome at a time.")
parser.add_argument("--output_dir", type=str, required=True, help="Directory to save entropy calculation results.")
parser.add_argument("--chromosome", type=str, default="chr1", help="Chromosome to analyze, default is chr1.")

args = parser.parse_args()

output_dir = args.output_dir
chromosomes = args.chromosome # can be space-separated list, or a single chromosome. 
                              # when multiple chormosomes are provided, the 1st chromosome will be used as the lead chromosome to load Tn5 insertion record, and the rest of the chromosomes will be used to update the Tn5 insertion record for each cell barcode.
chromosomes = chromosomes.split(' ') # turning into a list 
Nused_chromosomes = len(chromosomes)
#%%

# Step 1. Load Tn5 insertion frequency record for each cell barcode 
if Nused_chromosomes > 1:
    # When multiple chromosomes are used, 
    # combine the records for all the chromosomes by concatenating the Tn5 insertion frequency arrays for each cell barcode.
    lead_chrom = chromosomes[0]

    # Load Tn5 insertion record for each chromosomes     
    insert_record_dict = {}
    
    for chrom in chromosomes:
        tn5_insert_record_file = os.path.join(output_dir, f'{chrom}_insert_frequency.pickle')
        with open(tn5_insert_record_file, 'rb') as file:
            insert_record_dict[chrom] = pickle.load(file)
    
    # Prepare records for only the barcodes that have any fragments in the leading chromosome
    insert_record = defaultdict(return_none)
    
    keys = list(insert_record_dict[lead_chrom].keys())
    
    for bc in keys:
        chromosome_insert_arrays = [
            insert_record_dict[chrom].pop(bc) for chrom in chromosomes if insert_record_dict[chrom][bc] is not None
        ]
        if len(chromosome_insert_arrays) == Nused_chromosomes:
            insert_record[bc] = hstack(chromosome_insert_arrays) 
    
else:
    # When only one chromosome is used,
    # just load the Tn5 insertion record for that chromosome.
    chrom = chromosomes[0]
    
    tn5_insert_record_file = os.path.join(output_dir, f'{chrom}_insert_frequency.pickle')
    with open(tn5_insert_record_file, 'rb') as file:
        insert_record = pickle.load(file)



# %%
# Step 2. Calculate entropy for each cell barcode for the given chromosome 
barcode_entropy = {}

barcode_entropy_df_file = os.path.join(output_dir, f'calculated_barcode_entropy_df.tsv')


# Wrapper function for multiprocessing
def mp_wrapper(bc):
    Entropy_mixturedist, Entropy_open_region, Entropy_closeregion, Pi_closestate, mle_lambda, lambda_closeregion, p0, p_closed_state = mixdist_mle_entropy(insert_record[bc])
    return ( bc, Entropy_mixturedist, Entropy_open_region, Entropy_closeregion, Pi_closestate, mle_lambda, lambda_closeregion, p0, p_closed_state )

## << Move all the code before multiprocessing outside of __main__ 
## <<



# %%
if __name__ == "__main__":
    '''
    running example:
    python3 calculate_entropy.py --output_dir 'res' 
    
    # input arguments:
    # output_dir: dir to save entropy calculation results 
    # chromosome: chromosome to analyze, default is chr1    
    '''


    #%%
    ## Multiprocessing implementation
    if os.cpu_count() and os.cpu_count() >4:
        num_cpus = 4


        from multiprocessing import Pool       # start num_cpus worker processes
        keys = list(insert_record.keys())
        with Pool(processes=num_cpus)  as pool:
            res = pool.imap_unordered(mp_wrapper, keys, chunksize=4000)

            with open(barcode_entropy_df_file, 'w') as f:
                f.write(",Entropy,Entropy_open_region,Entropy_closeregion,Pi_closestate,mle_lambda,lambda_closeregion,P0_open_region,P_closed_state\n")
                for r in tqdm(res, total=len(keys), desc="Calculating entropy for each cell barcode"):
                    bc, Entropy_mixturedist, Entropy_open_region, Entropy_closeregion, Pi_closestate, mle_lambda, lambda_closeregion, p0, p_closed_state = r
                    if Entropy_mixturedist is None:
                        continue
                    barcode_entropy[bc] = Entropy_mixturedist
                    f.write(f"{bc},{Entropy_mixturedist},{Entropy_open_region},{Entropy_closeregion},{Pi_closestate},{mle_lambda},{lambda_closeregion},{p0},{p_closed_state}\n")


    ## Single-threaded implementation 
    else:
        with open(barcode_entropy_df_file, 'w') as f:
            f.write(",Entropy,Entropy_open_region,Entropy_closeregion,Pi_closestate,mle_lambda,lambda_closeregion,P0_open_region,P_closed_state\n")
            for bc, v in tqdm(insert_record.items(), desc="Calculating entropy for each cell barcode"):
                # Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = mixdist_mle_entropy(v)
                Entropy_mixturedist, Entropy_open_region, Entropy_closeregion, Pi_closestate, mle_lambda, lambda_closeregion, p0, p_closed_state = mixdist_mle_entropy(v)
                if Entropy_mixturedist is None:
                    continue
                barcode_entropy[bc] = Entropy_mixturedist
                f.write(f"{bc},{Entropy_mixturedist},{Entropy_open_region},{Entropy_closeregion},{Pi_closestate},{mle_lambda},{lambda_closeregion},{p0},{p_closed_state}\n")


    # save entropy results
    barcode_entropy_file = os.path.join(output_dir, f'calculated_barcode_entropy.pickle')
    with open(barcode_entropy_file, 'wb') as file:
        pickle.dump(barcode_entropy, file, protocol=pickle.HIGHEST_PROTOCOL)