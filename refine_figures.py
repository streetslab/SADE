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

    args = parser.parse_args()
    output_dir = args.output_dir
    chromosome = args.chromosome
    EntropyThreshold = args.entropythreshold
    
    import os 
    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')
    
    with open(entropy_file, 'rb') as f:
        entropy_data = pickle.load(f)
        
        
    entropy = np.array(list(entropy_data.values()))
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
    # %%