from multiprocessing import Pool
import time
import os 
import pickle
import numpy as np
from scipy.special import lambertw
import scipy.sparse
import pandas as pd
from tqdm import tqdm




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
    
    Entropy_mixturedist = Entropy_states + p_openstate * Entropy_openregion
    
    return Entropy_mixturedist, Entropy_openregion, Mnon0, mle_lambda, P0, p_closestate


def mp_wrapper(bc):

    Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = mixdist_mle_entropy(tn5_insert_freq[bc])
    return ( bc, Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state )


if __name__ == '__main__':
    ws = 3000

    output_dir = f'/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS{ws}'
    chromosome = 'chr1'


    tn5_insert_freq_file = os.path.join(output_dir, f'{chromosome}_insert_frequency.pickle')

    with open(tn5_insert_freq_file, 'rb') as f:
        tn5_insert_freq = pickle.load(f)
        
    
    results_df = pd.DataFrame( columns=[
        'entropy_mixdist_c',
        'entropy_open_region', 
        'Mp', 
        'mle_lambda', 
        'P0_open_region', 
        'P_closed_state'
    ] )

    entropy_df_file = os.path.join(output_dir, f'Multiprocessing_{chromosome}_barcode_entropy_df.csv') 
    with open(entropy_df_file, 'w') as f:
        
        f.write("Barcode,Entropy_mixturedist,Entropy_open_region,Mp,mle_lambda,P0_open_region,P_closed_state\n")
        
        with Pool(processes=4) as pool:         # start 4 worker processes
            keys = list(tn5_insert_freq.keys())
            res = pool.imap_unordered(mp_wrapper, keys, chunksize=4000)
            for r in tqdm(res, total=len(keys)):
                bc, Entropy_mixturedist, Entropy_open_region, Mnon0, mle_lambda, p0, p_closed_state = r
                if Entropy_mixturedist is None:
                    continue
                f.write(f"{bc},{Entropy_mixturedist},{Entropy_open_region},{Mnon0},{mle_lambda},{p0},{p_closed_state}\n")

        
        #print(results_df.head())
        
    # with Pool(processes=16) as pool:         # start 4 worker processes
    #     def f(x):
    #         return x*x
    #     result = pool.apply_async(f, (10,)) # evaluate "f(10)" asynchronously in a single process
    #     print(result.get(timeout=1))        # prints "100" unless your computer is *very* slow

    #     print(pool.map(f, range(10)))       # prints "[0, 1, 4,..., 81]"

    #     it = pool.imap_unordered(f, range(10), chunksize=2)
    #     print(next(it))                     # prints "0"
    #     print(next(it))                     # prints "1"
    #     print(next(it))                     # prints "4"
    #     print(next(it))                     # prints "9"
    #     print(it.next(timeout=1))           # prints "4" unless your computer is *very* slow


        
        
