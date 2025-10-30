#%% 
import sklearn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os



#%%
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--chromosome', type=str, default='chr1', help='Chromosome to process')
    
    
    
    args = parser.parse_args()
    output_dir = args.output_dir
    chromosome = args.chromosome
    
    #%%
    # Load entropy cutoff threshold
    entropy_cutoff_file = os.path.join(output_dir, 'entropy_cutoff.csv')
    with open(entropy_cutoff_file, 'r') as f:
        entropy_cutoff = float(f.readline().strip().split(',')[1])
    # entropy_cutoff = 10**(-1.5)
    # Load entropy data
    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')

    with open(entropy_file, 'rb') as f:
            entropy_data = pickle.load(f)

    entropy_df = pd.DataFrame.from_dict(entropy_data, orient='index', columns=['entropy'])
    
    # Load total fragments counts
    total_frag_counts_file = os.path.join(output_dir, 'total_fragments_counts.txt')
    total_frag_counts = pd.read_csv(total_frag_counts_file, sep=" ", index_col=1, header=None, names=['total_fragments'])

    
    #%% 
    all_bc_df = entropy_df.join(total_frag_counts, how='inner')  
    all_bc_df['log(total_fragments)'] = np.log10(all_bc_df['total_fragments'])
    all_bc_df['log(entropy)'] = np.log10(all_bc_df['entropy'])
    
    
    pass_entropy_df = all_bc_df.loc[all_bc_df['entropy'] >= entropy_cutoff].copy()
    
    


    #%% 
    # Linear regression
    from sklearn.linear_model import LinearRegression
    lr_fit = LinearRegression(fit_intercept=True).fit(X=pass_entropy_df[['log(total_fragments)']], y=pass_entropy_df['log(entropy)'])


    #%% 
    # Prediction for plot 
    min_log_total_fragments = pass_entropy_df['log(total_fragments)'].min()
    predict_log_entropy = lr_fit.predict(min_log_total_fragments.reshape(-1, 1))[0]
    
    #%% 
    # Plot
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(all_bc_df['log(total_fragments)'], all_bc_df['log(entropy)'], s=1, color='gray', label='All barcodes')
    ax.set_xlabel('log10(Total Fragments)')
    ax.set_ylabel('log10(Entropy)')
    ax.scatter(pass_entropy_df['log(total_fragments)'], pass_entropy_df['log(entropy)'], s=1, color='blue', label='Pass Entropy')
    ax.set_xlabel('log2(Total Fragments)')
    ax.axline([min_log_total_fragments, predict_log_entropy], slope=lr_fit.coef_[0], color='red', linestyle='--', label='Linear Regression')
    plt.tight_layout()
    plt.show()
    
# %%
output_dir = '/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000'
# output_dir = '/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS30000'
# output_dir = '/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS300'

# To select range of total fragments to explot on WS's impact on their entropy values differences between good and bad barcodes
min_log_total_fragments = 2.7
max_log_total_fragments = 3


# %%
