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


    entropy_all_file = os.path.join(output_dir, f'calculated_barcode_entropy_df.tsv')
    entropy_all_df = pd.read_csv(entropy_all_file, sep=',', header=0, index_col=0)
    
    entropy_w_DNAdebrisflag_file = os.path.join(output_dir, 'entropy_filtered_bc_w_DNAdebrisflag_df.tsv')    
    entropy_w_DNAdebrisflag_df = pd.read_csv(entropy_w_DNAdebrisflag_file, sep=',', index_col=0)

    entropy_df_file = os.path.join(output_dir, f'entropy_filtered_bc_df.tsv')
    entropy_df = pd.read_csv(entropy_df_file, sep=',', header=0, index_col=0)
 

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
    frag_df.columns = ['chrom', 'start', 'end', 'atac_barcodes', 'support']
    frag_df['length'] = frag_df['end'] - frag_df['start']

    frag_df['entropy_pass'] = frag_df['atac_barcodes'].isin(entropy_df.index)
    filter_frag_df = frag_df[frag_df['entropy_pass']]
    garbage_frag_df = frag_df[~ frag_df['entropy_pass']]

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 6), sharey=True)
    ax[0].hist(filter_frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[0].set_xlabel('Fragment length')
    ax[0].set_ylabel('Frequency')
    ax[0].set_title('Reads in permitted droplets')
    ax[1].hist(frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[1].set_xlabel('Fragment length')
    ax[1].set_ylabel('Frequency')
    ax[1].set_title('All reads in the experiment')
    ax[2].hist(garbage_frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[2].set_xlabel('Fragment length')
    ax[2].set_ylabel('Frequency')
    ax[2].set_title('Reads in filtered out droplets')
    fig_file = os.path.join(output_dir, 'figures', 'frag_length_distribution.png')
    fig.tight_layout()
    fig.savefig(fig_file)

    #%%
    # Demonstration figure 
    total_fragments = frag_df.groupby('atac_barcodes').size().reset_index(name='total_fragments')
    total_fragments.set_index('atac_barcodes', inplace=True)
    
    for df in [entropy_all_df, entropy_w_DNAdebrisflag_df, entropy_df]:
        df['log_total_fragments'] = np.log1p(df.index.map(total_fragments['total_fragments']))
        df['log_entropy'] = np.log10(df['Entropy'])
    
    #%%
    fig = plt.figure(layout="constrained", figsize=(6.6, 5))

    # Create a top-level GridSpec with 3 rows
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1])

    # --- Row 1 ---
    gs0 = gs[0].subgridspec(1, 2, width_ratios=[3, 8])
    ax1 = fig.add_subplot(gs0[0, 0])
    ax2 = fig.add_subplot(gs0[0, 1])

    # --- Row 2 ---
    gs1 = gs[1].subgridspec(1, 2, width_ratios=[3, 8])
    ax3 = fig.add_subplot(gs1[0, 0]) 
    ax4 = fig.add_subplot(gs1[0, 1])

    # --- Row 3 ---
    # ax5 = fig.add_subplot(gs[2, 0])
    gs2 = gs[2].subgridspec(1, 2, width_ratios=[3, 8])
    ax5 = fig.add_subplot(gs2[0, 0])
    ax6 = fig.add_subplot(gs2[0, 1])


    entropy_min, entropy_max = entropy_all_df['log_entropy'].min(), entropy_all_df['log_entropy'].max()

    # ax1
    s1 = ax1.hist(entropy_all_df['log_entropy'], bins=140, alpha=0.5, color='black', 
                orientation='horizontal')
    ax1.set_ylim(entropy_min - 0.5   , entropy_max + 0.5)
    ax1.invert_xaxis()
    # ax1.yaxis.tick_right()
    ax1.set_xlabel('Frequency')
    ax1.set_ylabel('Entropy ($log$)')
    # ax2 
    sorted_log10entropy_list = entropy_all_df['log_entropy'].dropna().sort_values(ascending=False).values
    rank_list = np.arange( len(sorted_log10entropy_list) )
    entropy_cutoff =   entropy_df['log_entropy'].min()
    s2 = ax2.plot(rank_list, sorted_log10entropy_list,  'o', markersize=3, color='grey')
    ax2.axhline(entropy_cutoff, color='blue', linestyle='--', alpha=0.7)
    ax2.set_xlabel('Rank')
    # set ax2 x-ticks to be none
    ax2.set_yticks([])

    # ax3 
    DNAdebris_index = entropy_w_DNAdebrisflag_df['DNA_debris'] == 'YES'
    noDNAdebris_index = entropy_w_DNAdebrisflag_df['DNA_debris'] == 'NO'
    s3 = ax3.hist(entropy_w_DNAdebrisflag_df['log_entropy'][noDNAdebris_index], bins=50, alpha=0.6, color="#2E5BFF",
                orientation='horizontal')
    s3_a = ax3.hist(entropy_w_DNAdebrisflag_df['log_entropy'][DNAdebris_index],  alpha=0.6, color="#E4D00A",
                orientation='horizontal')
    ax3.set_ylim(entropy_cutoff - 0.5, entropy_max + 0.5)
    ax3.invert_xaxis()
    ax3.set_xlim(ax1.get_xlim()[::-1]) # make x-axis of ax3 same as ax1
    ax3.set_ylabel('Entropy ($log$)')
    ax3.set_xlabel('Frequency')

    # ax4 
    s4 = ax4.scatter(entropy_w_DNAdebrisflag_df['log_total_fragments'], entropy_w_DNAdebrisflag_df['log_entropy'],
                    c = 1- entropy_w_DNAdebrisflag_df['P_closed_state'],
                    alpha=1, s=5, cmap='coolwarm' )
    ax4.set_xlabel('Total Fragments ($log$)')
    ax4.set_yticklabels([]) # hide y-ticks
    ax4.set_ylim(entropy_cutoff - 0.5, entropy_max + 0.5)


    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # Create colorbar axes OUTSIDE ax2 [x_offset, y_offset, width, height]
    cax2 = inset_axes(ax4,
                    width="3%",  
                    height="100%",
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0., 1, 1), # (x, y, width, height) relative to ax2
                    bbox_transform=ax4.transAxes,
                    borderpad=0)
    fig.colorbar(s4, cax=cax2, label='$P_{open-state}$')
    # ax5 
    ax5.axis(False)
    # ax6
    s6 = ax6.scatter(entropy_all_df['log_total_fragments'], entropy_all_df['log_entropy'], 
                    c ='grey', alpha=0.4, s=3)
    ax6.set_xlabel('Total Fragments ($log$)')
    ax6.set_ylabel('Entropy ($log$)')
    ax6.scatter(entropy_df['log_total_fragments'], entropy_df['log_entropy'], 
                    c ='#2E5BFF', alpha=0.4, s=3)


# %%
