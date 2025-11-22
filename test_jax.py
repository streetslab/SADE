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

import dask.array as da
import dask
import sparse

precision = da.array(1e-7) # avoid nan in lambert W function calculation

# %%

# For loading pickle files of Tn5 insertion frequency that has defaultdict structure
from utils import return_none


# %%

    
output_dir = '/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000_mixD'
chromosome = 'chr1'

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

#%%
da_record = {}
i = 0
for bc, v in tqdm(insert_record.items(), desc="Calculating entropy for each cell barcode"):
    i += 1
    # Entropy_mixturedist, Entropy_open_region, Entropy_open_region_noadj, mle_lambda, P0, p_closed_state, divergence = mixdist_mle_entropy(v)
    v_r, v_c, v_data = scipy.sparse.find(v)
    d_v = da.from_array(v_data, chunks='auto')
    da_record[bc] = d_v
    
    # if i == 100:
    #     break
    # barcode_entropy[bc] = Entropy_mixturedist
    # entropy_df.loc[bc, 'Entropy'] = Entropy_mixturedist
    # entropy_df.loc[bc, 'Entropy_open_region'] = Entropy_open_region
    # entropy_df.loc[bc, 'Entropy_open_region_noadj'] = Entropy_open_region_noadj
    # entropy_df.loc[bc, 'mle_lambda'] = mle_lambda
    # entropy_df.loc[bc, 'P0_open_region'] = P0
    # entropy_df.loc[bc, 'P_closed_state'] = p_closed_state
    # entropy_df.loc[bc, 'kl_divergence'] = divergence
# %%
# v_r, v_c, v_data = scipy.sparse.find(v)
# d_v = da.from_array(v_data)

def mixdist_mle_entropy(d_v:da.array, Mtotal:da.array):

    M_non0 = da.array(d_v.shape[0])

    adjust_v = (d_v + 1) // 2 # adjust for paired-end sequencing
    cur_adj, freq_adj = da.unique(adjust_v, return_counts=True)

    ZTP_prob = freq_adj / M_non0 # zero-truncated Poisson probability for each unique insertion count
    mean_adj = da.matmul( cur_adj, ZTP_prob ) 

    # Lambert W function requires (-M1 * exp(-M1)) >= -1/e  --> mean_adj > 1
    lambda_mle = lambertw( -mean_adj * da.exp(-mean_adj) , k=0).real + mean_adj
    P0 = da.exp(-lambda_mle) 
    
    N0_openregion = M_non0 / ( 1- P0) * P0 # integer restriction does not matter in calculation of entropy
    
    # Entropy of mix states
    p_closestate = N0_openregion / Mtotal
    Entropy_states = -da.log2( p_closestate ) * p_closestate - da.log2(1 - p_closestate) * (1 - p_closestate) #TODO:  make it vector multiplication
    
    # Entropy only for open regions under Poisson model
    Poisson_prob = ZTP_prob * (1 - P0) # Only for non-zero windows in open region
    Entropy_openregion = -da.matmul(Poisson_prob, da.log2(Poisson_prob)) - P0 * da.log2(P0)  # including imputed zero-insertion windows in open region

    Entropy_mixturedist = Entropy_states + (1 - p_closestate) * Entropy_openregion

    return Entropy_mixturedist, Entropy_openregion, P0, p_closestate

# %%


def calcualte_P0(d_v:da.array):

    M_non0 = da.array(d_v.shape[0])

    adjust_v = (d_v + 1) // 2 # adjust for paired-end sequencing
    cur_adj, freq_adj = da.unique(adjust_v, return_counts=True)

    ZTP_prob = freq_adj / M_non0 # zero-truncated Poisson probability for each unique insertion count
    mean_adj = da.matmul( cur_adj, ZTP_prob ) 

    # Lambert W function requires (-M1 * exp(-M1)) >= -1/e  --> mean_adj > 1
    lambda_mle = lambertw( -mean_adj * da.exp(-mean_adj) , k=0).real + mean_adj
    P0 = da.exp(-lambda_mle) 
    
    return P0


#%% 
P0_result = []
for k, v in tqdm(da_record.items(), desc="Calculating P0 for each cell barcode with Dask"):
    P0 = calcualte_P0(v)
    P0_result.append(P0)
    
P0_result_computed = dask.compute( *P0_result )


#%%
from dask.diagnostics import ProgressBar
Mtotal = insert_record[bc].shape[1]
entropy_mixdist = []
with ProgressBar():
    for k, v in tqdm(da_record.items(), desc="Calculating entropy for each cell barcode with Dask"):
        Entropy_mixturedist, Entropy_open_region, P0, p_closed_state = mixdist_mle_entropy(v, Mtotal=da.array(Mtotal))
        entropy_mixdist.append( Entropy_mixturedist )
        

# %%
# from dask.distributed import Client, progress
# client = Client(threads_per_worker=4, n_workers=1)
# client
#%%
pbar = ProgressBar()
pbar.register()

final_results = dask.compute( *entropy_mixdist )

pbar.unregister()
# %%
