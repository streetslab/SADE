#%%
#### GLOBAL VARIABLES >>
precision = 1e-10 # for plot
#### GLOBAL VARIABLES << 

# %%

import pandas as pd
import os
import numpy as np

import matplotlib.pyplot as plt

from utils import get_frag_overlap_peaks_df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot entropy distribution for cell barcodes and filtered fragments.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument('--chromosome', type=str, default='chr1', help='Chromosome used to calculate entropy.')
    parser.add_argument("--crbarcode_file", type=str, required=True, help="CellRanger barcode file, i.e., barcodes labeled to contain cells.")

    
    args = parser.parse_args()
    output_dir = args.output_dir
    crbarcode_file = args.crbarcode_file
    chromosome = args.chromosome

    # %%

    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy_df.tsv')
    entropy_df = pd.read_csv(entropy_file, sep='\t', header=0, index_col=0)
    
    entropycutoff_file = os.path.join(output_dir, 'entropy_cutoff.csv')
    with open(entropycutoff_file, 'r') as f:
        line = f.readline()
        entropythreshold = float(line.strip().split(',')[1])
        
    entropy_df['log10_entropy'] = np.log10(entropy_df['Entropy'] + precision)
    entropy_df['pass'] = entropy_df['Entropy'] > entropythreshold
    

    figure_subdir = os.path.join(output_dir, 'figures')
    if not os.path.exists(figure_subdir):
        os.makedirs(figure_subdir)
    # %%

    # >>>  compare entropy distribution between CR-labeled cells vs empty barcodes <<<
    cr_bc = pd.read_csv(crbarcode_file, sep='\t', header=None)

    entropy_df['cr_cell'] = entropy_df.index.isin(cr_bc[0])


    cr_entropy_df = entropy_df.loc[entropy_df['cr_cell']]
    cr_empty_entropy = entropy_df[~entropy_df['cr_cell']]

    print("CellRanger labeled good vs bad barcodes stats:")
    print(f"    Number of good barcodes: {len(cr_entropy_df)}")
    print(f"    Number of bad barcodes: {len(cr_empty_entropy)}")
    print(f"    Number of barcodes : {len(entropy_df)} \n")

    with open(os.path.join(output_dir, 'entropy_stats.txt'), 'w') as f:
        f.write("Entropy statistics for CR_labeled cells:\n")
        cr_entropy_df['Entropy'].describe().to_string(f)
        f.write('\n')
        f.write('entropythreshold: ' + str(entropythreshold) + '\n')
        f.write('\n')
        f.write("Entropy statistics for CR_empty barcodes:\n")
        cr_empty_entropy['Entropy'].describe().to_string(f)
        
    print("Entropy filtering good vs bad barcodes stats: ")
    print(f"    Number of good barcodes : {entropy_df['pass'].sum()}")
    print(f"    Number of bad barcodes : {len(entropy_df) - entropy_df['pass'].sum()}")
    print(f"    Total number of barcodes : {len(entropy_df)} \n")


    # %%
    # remove NaN values for plot
    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].violinplot(cr_entropy_df['Entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[1].violinplot(cr_empty_entropy['Entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[0].set_ylabel('Entropy')
    ax[0].set_xlabel('CR cell barcodes')
    ax[1].set_xlabel('CR empty barcodes')
    ax[0].axhline(y=entropythreshold, color='r', linestyle='-', label=f'{entropythreshold:.4f}')
    ax[1].axhline(y=entropythreshold, color='r', linestyle='-', label=f'{entropythreshold:.4f}')
    fig.legend()
    fig_file = os.path.join(figure_subdir, 'entropy_violin_EmptyCRcell.png')
    fig.savefig(fig_file)
    
    
    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].violinplot(cr_entropy_df['log10_entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[1].violinplot(cr_empty_entropy['log10_entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[0].set_ylabel('log10 Entropy')
    ax[0].set_xlabel('CR cell barcodes')
    ax[1].set_xlabel('CR empty barcodes')
    ax[0].axhline(y=np.log10(entropythreshold ), color='r', linestyle='-', label=f'log10({entropythreshold:.4f})')
    ax[1].axhline(y=np.log10(entropythreshold ), color='r', linestyle='-', label=f'log10({entropythreshold:.4f})')
    fig.legend()
    fig_file = os.path.join(figure_subdir, 'log10_entropy_violin_EmptyvsCell.png')
    fig.savefig(fig_file)


    #%%
    
    entropy_df_sort = entropy_df.sort_values(by='Entropy', ascending=False)
    entropy_df_sort['index'] = np.arange(1, entropy_df_sort.shape[0] + 1)
    
    # Knee plot for entropy and colored by CellRanger cell-calling
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(entropy_df_sort['index'], np.log10(entropy_df_sort['Entropy'] + precision), \
                c=entropy_df_sort['cr_cell'].map({True: 'blue', False: 'orange'}), \
                marker= '.', alpha=0.4, label='CR Cell Barcodes', s=3)
    ax.axhline(y=np.log10(entropythreshold + precision), color='green', linestyle='--', label=f'Entropy cutoff: log10({entropythreshold:.4f})', lw=2)
    ax.set_ylabel('Entropy (log10-scaled)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_color_crbc.png')
    fig.savefig(fig_file, bbox_inches='tight')

    #%%
    Nbcs_to_plot = 30000  # 10x can only process up to 20k cells
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(entropy_df_sort['index'][:Nbcs_to_plot], np.log10(entropy_df_sort['Entropy'][:Nbcs_to_plot] + precision), \
                c=entropy_df_sort['cr_cell'][:Nbcs_to_plot].map({True: 'blue', False: 'orange'}), \
                marker= '.', alpha=0.4, label='CR Cell Barcodes', s=1)
    ax.axhline(y=np.log10(entropythreshold + precision), color='green', linestyle='--', label=f'log10({entropythreshold:.4f})', lw=2)
    ax.set_ylabel('Entropy (log10-scaled)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_color_crbc_top30k.png')
    fig.savefig(fig_file, bbox_inches='tight')




    # <<< Compare <<<


