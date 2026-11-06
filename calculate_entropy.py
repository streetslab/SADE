#%% 

import os
import re
import sys 
import numpy as np

import scipy

import pickle

from tqdm import tqdm


from scipy.special import lambertw
from scipy.stats import poisson 
from scipy.sparse import hstack, lil_matrix
from collections import defaultdict

import argparse

#%%

from utils import return_none, get_chromsize_dict, default_chrome_bin_sparse



#%%
############################################# 
## Calculate entropy for each cell barcode >>
#############################################
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
    # Entropy_mixturedist = Entropy_states + p_openstate * Entropy_openregion / (Mnon0 * mean_adj) **2  # Adjust for sample size effect
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
parser.add_argument("--genomesize_file", type=str, required=True, help="Genome size reference file.")
parser.add_argument("--shortcut_mode", type=str, required=False,  default="no", help="Go straight to Tn5 insertion counting and skip fragment splitting, default is no. no or yes")
parser.add_argument("--nprocess", type=str, required=False, help="Number of cpu cores to use.")


args = parser.parse_args()

output_dir = args.output_dir
chromosomes_spec = args.chromosome # can be space-separated list, or a single chromosome. 
                              # when multiple chormosomes are provided, the 1st chromosome will be used as the lead chromosome to load Tn5 insertion record, and the rest of the chromosomes will be used to update the Tn5 insertion record for each cell barcode.
# chromosomes = chromosomes_spec.split(' ') # turning into a list 
genome_chromsize_file = args.genomesize_file
chromsize_dict = get_chromsize_dict(genome_chromsize_file)


try:
    nprocess = int(args.nprocess) if args.nprocess is not None else 4
except ValueError:
    nprocess = 4 # default number of CPU cores to use if the provided value is invalid

shortcut_mode = args.shortcut_mode.lower() == "yes"
if not shortcut_mode:
    save_subdir = output_dir
else:
    save_subdir = os.path.join(output_dir, 'per_chromosome_entropy')
    if not os.path.exists(save_subdir):
        os.makedirs(save_subdir, exist_ok=True)


#%%
# Step 1. Load Tn5 insertion frequency record for each cell barcode 
if len(chromosomes_spec.split(' ')) == 1 and chromosomes_spec != "autosomes" :
    # When only one chromosome is used,
    # just load the Tn5 insertion record for that chromosome.
    chromosomes = chromosomes_spec.split(' ') # turning into a list
    chrom = chromosomes[0]
    
    tn5_insert_record_file = os.path.join(save_subdir, f'{chrom}_insert_frequency.pickle')
    with open(tn5_insert_record_file, 'rb') as file:
        insert_record = pickle.load(file)
else:
    # When multiple chromosomes are used, 
    # combine the records for all the chromosomes by concatenating the Tn5 insertion frequency arrays for each cell barcode.
    if chromosomes_spec.lower() == "autosomes":
        # all autosomes used 
        chromosomes = []
        for chrom in chromsize_dict.keys():
            if re.match(r"^chr[0-9]+$", chrom) is not None:
                chromosomes.append(chrom)
                
    else:
        chromosomes = chromosomes_spec.split(' ')
    
    lead_chrom = chromosomes[0]
    Nchromosomes = len(chromosomes)
    
    # Load Tn5 insertion record for each chromosomes     
    insert_record_dict = {}
    keys = set()
    for chrom in chromosomes:
        tn5_insert_record_file = os.path.join(save_subdir, f'{chrom}_insert_frequency.pickle')
        with open(tn5_insert_record_file, 'rb') as file:
            insert_record_dict[chrom] = pickle.load(file)
            keys.update(list(insert_record_dict[chrom].keys()))
    
    # Prepare records for only the barcodes that have any fragments in the leading chromosome
    insert_record = defaultdict(return_none)
    
    chrom_bin_size  = {}
    for chrom in chromosomes:
        random_bc = list(insert_record_dict[chrom].keys())[0]
        chrom_bin_size[chrom] = insert_record_dict[chrom][random_bc].shape[1]

    for bc in keys:
        chromosome_insert_arrays = [
            insert_record_dict[chrom].pop(bc) if insert_record_dict[chrom][bc] is not None else lil_matrix((1, chrom_bin_size[chrom]), dtype=np.int32)  for chrom in chromosomes
        ] # fill in empty insertions when a barcode has no fragments in a particular chromosome
        
        insert_record[bc] = hstack(chromosome_insert_arrays) 



# %%
# Step 2. Calculate entropy for each cell barcode for the given chromosome 
barcode_entropy = {}
if not shortcut_mode:
    barcode_entropy_df_file = os.path.join(save_subdir, f'calculated_barcode_entropy_df.tsv')
else:
    # Always specific chromosome is specified in shortcut mode
    if chromosomes_spec.lower() == "autosomes" or len(chromosomes) > 1:
        sys.exit("Shortcut mode with 'autosomes' or multiple chromosomes is not supported in per-chromosome calculation.")
    else:
        barcode_entropy_df_file = os.path.join(save_subdir, f'{chromosomes[0]}_calculated_barcode_entropy_df.tsv')

# Wrapper function for multiprocessing
def mp_wrapper(bc):
    Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = mixdist_mle_entropy(insert_record[bc])
    return ( bc, Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state )

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
    if os.cpu_count() and os.cpu_count() > nprocess:
        num_cpus = nprocess


        from multiprocessing import Pool       # start num_cpus worker processes
        keys = list(insert_record.keys())
        with Pool(processes=num_cpus)  as pool:
            res = pool.imap_unordered(mp_wrapper, keys, chunksize=4000)

            with open(barcode_entropy_df_file, 'w') as f:
                f.write(",Entropy,Entropy_open_region,Mp,mle_lambda,P0_open_region,P_closed_state\n")
                for r in tqdm(res, total=len(keys), desc="Calculating entropy for each cell barcode"):
                    bc, Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = r
                    if Entropy_mixturedist is None:
                        continue
                    barcode_entropy[bc] = Entropy_mixturedist
                    f.write(f"{bc},{Entropy_mixturedist},{Entropy_open_region},{Mnon0},{mle_lambda},{p0},{p_closed_state}\n")


    ## Single-threaded implementation 
    else:
        with open(barcode_entropy_df_file, 'w') as f:
            f.write(",Entropy,Entropy_open_region,Mp,mle_lambda,P0_open_region,P_closed_state\n")
            for bc, v in tqdm(insert_record.items(), desc="Calculating entropy for each cell barcode"):
                Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = mixdist_mle_entropy(v)
                if Entropy_mixturedist is None:
                    continue
                barcode_entropy[bc] = Entropy_mixturedist
                f.write(f"{bc},{Entropy_mixturedist},{Entropy_open_region},{Mnon0},{mle_lambda},{p0},{p_closed_state}\n")


    # save entropy results
    barcode_entropy_file = os.path.join(save_subdir, f'calculated_barcode_entropy.pickle')
    with open(barcode_entropy_file, 'wb') as file:
        pickle.dump(barcode_entropy, file, protocol=pickle.HIGHEST_PROTOCOL)