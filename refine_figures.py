

#%%


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import argparse
import os 


#%%


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Refine figures.")
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for figures.')
    parser.add_argument('--chromosome', type=str, default='chr1', help='Chromosome to analyze.')


    args = parser.parse_args()
    output_dir = args.output_dir
    chromosome = args.chromosome

    #%%


    figure_subdir = os.path.join(output_dir, 'figures')

    fragment_file = os.path.join(output_dir, 'fragments.tsv')

    entropy_cutoff_file = os.path.join(output_dir, 'entropy_cutoff.csv')
    with open(entropy_cutoff_file, 'r') as f:
        line = f.readline()
        EntropyThreshold = float(line.strip().split(',')[1])


    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy_df.tsv')
    entropy_df = pd.read_csv(entropy_file, sep=',', header=0, index_col=0)
    entropy_df['pass'] = entropy_df['Entropy'] > EntropyThreshold

    entropy = entropy_df['Entropy'].dropna().values # remove None values   
    entropy.sort() # sort the entropy values


    # %%
    # Knee plot for entropy
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(np.log10(entropy[::-1] ), '-', color='blue', alpha=0.5, label='Entropy')
    ax.axhline(y=np.log10(EntropyThreshold ), color='green', linestyle='--', label=f'log10({EntropyThreshold:.4f})')
    ax.set_ylabel('Entropy (log10)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_logscale.png')
    fig.savefig(fig_file, bbox_inches='tight')
    #%%
    #histogram of log10(entropy) values -- for (potentially) Gaussian mixture model fitting
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(np.log10(entropy ), bins=100, color='blue',  label='log10(Entropy)')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('log10(Entropy)')
    # fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_Entropy_histogram.png')
    fig.savefig(fig_file, bbox_inches='tight')

    # histogram of log10(OpenRegion-entropy) values 
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(np.log10(entropy_df['Entropy_open_region']/ entropy_df['Mp']), bins=100, color='blue',  log=False, label='log10(Open-region Entropy)')
    ax.set_ylabel('log10(Frequency)')
    ax.set_xlabel('log10(Open-region Entropy)')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_OpenregionEntropy_histogram.png')
    fig.savefig(fig_file, bbox_inches='tight')

    # %%
    # Distribution of fragments length 
    frag_df = pd.read_csv(fragment_file, sep='\t', header=None)
    frag_df.columns = ['chrom', 'left_insert', 'right_insert', 'cb', 'support']
    frag_df['length'] = frag_df['right_insert'] - frag_df['left_insert']

    frag_df['entropy_pass'] = frag_df['cb'].isin(entropy_df[entropy_df['pass']].index)
    filter_frag_df = frag_df[frag_df['entropy_pass']]
    garbage_frag_df = frag_df[~ frag_df['entropy_pass']]

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 6), sharey=True)
    ax[0].hist(filter_frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[0].set_xlabel('Fragment length')
    ax[0].set_ylabel('Frequency')
    ax[0].set_title('Fragment length distribution for filtered fragments')
    ax[1].hist(frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[1].set_xlabel('Fragment length')
    ax[1].set_ylabel('Frequency')
    ax[1].set_title('Fragment length distribution for all fragments')
    ax[2].hist(garbage_frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[2].set_xlabel('Fragment length')
    ax[2].set_ylabel('Frequency')
    ax[2].set_title('thrown-away fragments')
    fig_file = os.path.join(output_dir, 'figures', 'frag_length_distribution.png')
    fig.tight_layout()
    fig.savefig(fig_file)

    #%%
