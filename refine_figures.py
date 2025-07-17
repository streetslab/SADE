#%%
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%%

precision = 1e-10  # to avoid log(entropy==0) issues



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Refine figures.")
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for figures.')
    parser.add_argument('--chromosome', type=str, default='chr1', help='Chromosome to analyze.')
    parser.add_argument('--entropythreshold', type=float, default=0.001, help='Entropy threshold for filtering.')
    parser.add_argument("--crbarcode_file", type=str, required=True, help="CellRanger barcode file, i.e., barcodes labeled to contain cells.")


    args = parser.parse_args()
    output_dir = args.output_dir
    chromosome = args.chromosome
    EntropyThreshold = args.entropythreshold
    crbarcode_file = args.crbarcode_file
    
    import os 
    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')
    
    with open(entropy_file, 'rb') as f:
        entropy_dict = pickle.load(f)
        
        
    entropy = np.array(list(entropy_dict.values()))
    entropy.sort() # sort the entropy values

    # %%
    # Knee plot for entropy
    figure_subdir = os.path.join(output_dir, 'figures')
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.log10(entropy[::-1] + precision), '-', color='blue', alpha=0.5, label='Entropy')
    ax.axhline(y=np.log10(EntropyThreshold + precision), color='red', linestyle='--', label=f'log10("{EntropyThreshold}")')
    ax.set_ylabel('Entropy (log10-scaled)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_logscale.png')
    fig.savefig(fig_file, bbox_inches='tight')
    
    
    #%%
    #histogram of log10(entropy) values -- for (potentially) Gaussian mixture model fitting
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(np.log10(entropy + precision), bins=100, color='blue',  label='log10(Entropy)')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('log10(Entropy)')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_histogram.png')
    fig.savefig(fig_file, bbox_inches='tight')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(np.log10(entropy + precision), bins=100, color='blue',  log=True, label='log10(Entropy)')
    ax.set_ylabel('log10(Frequency)')
    ax.set_xlabel('log10(Entropy)')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_histogram_logFreq.png')
    fig.savefig(fig_file, bbox_inches='tight')
    # %%
    
    bc_entropy = pd.DataFrame.from_dict(entropy_dict, orient='index', columns=['entropy'])
    
    cr_bc = pd.read_csv(crbarcode_file, sep='\t', header=None)
    cr_bc_pd = cr_bc.set_index(0)
    cr_bc_pd['cr_cell'] = True

    bc_entropy_sort = bc_entropy.sort_values(by='entropy', ascending=False)
    bc_entropy_sort['index'] = np.arange(bc_entropy_sort.shape[0])
    bc_entropy_sort = bc_entropy_sort.join(cr_bc_pd, how='left')
    bc_entropy_sort['cr_cell'] = bc_entropy_sort['cr_cell'].fillna(False)
    
    # Knee plot for entropy and colored by CellRanger cell-calling 
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(bc_entropy_sort['index'], np.log10(bc_entropy_sort['entropy'] + precision), \
                c=bc_entropy_sort['cr_cell'].map({True: 'blue', False: 'orange'}), \
                marker= '.', alpha=0.4, label='CR Cell Barcodes', s=1)
    ax.axhline(y=np.log10(EntropyThreshold + precision), color='red', linestyle='--', label=f'log10({EntropyThreshold})')
    ax.set_ylabel('Entropy (log10-scaled)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_color_crbc.png')
    fig.savefig(fig_file, bbox_inches='tight')