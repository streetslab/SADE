#%%
import pandas as pd
import numpy as np
import os

#%%
MininumFragments = 1000

from utils import get_frag_overlap_peaks_df


#%%
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Visualize fragment overlaps with peaks.")
    parser.add_argument("--res_dir", required=True, type=str, help="Directory to results.")
    parser.add_argument('--chromosome', type=str, default='chr1', help='Chromosome to analyze.')
    
    args = parser.parse_args()
    output_dir = args.res_dir
    chromosome = args.chromosome
    
    #%%
    frag_overlap_peaks_df, frag_overlap_entropypeaks_df =  get_frag_overlap_peaks_df(output_dir, chromosome)
    # Save the processed dataframes as summary statistics 
    fragments_overlap_subdir = os.path.join(output_dir, 'fragments_overlap_peaks')
    
    frag_overlap_peaks_df.to_csv(os.path.join(fragments_overlap_subdir, 'summary_statistics_frag_overlap_peaks.csv'), sep='\t')
    frag_overlap_entropypeaks_df.to_csv(os.path.join(fragments_overlap_subdir, 'summary_statistics_frag_overlap_entropypeaks.csv'), sep='\t')


    #%%
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D


    # Filter out barcodes with less than MininumFragments fragments
    overlap_peaks_df = frag_overlap_peaks_df[frag_overlap_peaks_df['total_fragments'] > MininumFragments].sort_values(by='frag_overlap_peaks%', ascending=False)
    overlap_entropy_peaks_df = frag_overlap_entropypeaks_df[frag_overlap_entropypeaks_df['total_fragments'] > MininumFragments].sort_values(by='frag_overlap_entropy_peaks%', ascending=False)


    figure_subdir = os.path.join(output_dir, 'figures')
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.sort(overlap_peaks_df['frag_overlap_peaks%'].values)[::-1], label='Overlap with Peaks', color='orange', marker='.', markersize=2, alpha=0.5)
    ax.plot(np.sort(overlap_entropy_peaks_df['frag_overlap_entropy_peaks%'].values)[::-1], label='Overlap with Entropy Peaks', color='blue', marker='.', markersize=2, alpha=0.5)
    ax.set_xlabel('Rank of barcodes ')
    ax.set_ylabel('Percentage of fragments')
    fig.legend()
    fig_file = os.path.join(figure_subdir, 'fragments_overlap_peaks.png')
    fig.savefig(fig_file, dpi=800, bbox_inches='tight')

    #%% 
    fig, ax = plt.subplots(figsize=(10, 6))
    color_map = {'Overlap with Peaks':'orange', 'Overlap with Entropy Peaks':'blue'}
    ax.scatter(
        np.log10(frag_overlap_peaks_df['total_fragments'].values),
        frag_overlap_peaks_df['frag_overlap_peaks%'].values,
        s=2, alpha=0.5, label='Overlap with Peaks', color='orange'
    )
    ax.scatter(
        np.log10(frag_overlap_entropypeaks_df['total_fragments'].values),
        frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'].values,
        s=2, alpha=0.5, label='Overlap with Entropy Peaks', color='blue'
    )
    ax.set_xlabel('Log10 of total fragments', fontsize=14)
    ax.set_ylabel(f'Fragments% overlapping with peaks', fontsize=14)
    fig.legend(handles=[Line2D([0], [0], marker='o', color='w', label='Overlap with Peaks',
                          markerfacecolor='orange', markersize=10),
                 Line2D([0], [0], marker='o', color='w', label='Overlap with Entropy Peaks',
                           markerfacecolor='blue', markersize=10)],
           bbox_to_anchor=(1.02, 1), loc='upper left')
    fig_file = os.path.join(figure_subdir, 'fragments_overlap_peaks_scatter.png')
    fig.savefig(fig_file, dpi=800, bbox_inches='tight')


# %%
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(
        np.log10(frag_overlap_entropypeaks_df['total_fragments'].values),
        frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'].values,
        s=4, 
        c=np.log10(frag_overlap_entropypeaks_df['entropy'].values),
        cmap='viridis',
    )
    ax.set_xlabel('log10(Total Fragments)')
    ax.set_ylabel('Fragment Overlap Entropy Peaks %')
    ax.set_title('Fragments Overlap with Entropy Peaks Colored by Entropy')
    fig.colorbar(ax.collections[0], ax=ax, label='Log10 Entropy')
    fig_file = os.path.join(figure_subdir, 'fragments_overlap_entropypeaks_colorentropy.png')
    fig.savefig(fig_file, dpi=800, bbox_inches='tight')
# %%

    fig, ax = plt.subplots(figsize=(10, 6))
    p = ax.scatter(np.log10(frag_overlap_entropypeaks_df['total_fragments'].values), 
                   np.log10(frag_overlap_entropypeaks_df['entropy'].values), 
                   c=frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'], \
        cmap='viridis', alpha=0.6, s=10)
    fig.colorbar(p, ax=ax, label='Fragment Overlap Entropy Peaks %')
    ax.set_xlabel('log10(Total Fragments)')
    ax.set_ylabel('log10(Entropy)')
    fig_file = os.path.join(figure_subdir, f'log10_entropy_vs_total_fragments.png')
    fig.savefig(fig_file)

