#%%
import pandas as pd
import numpy as np
import os

#%%
MininumFragments = 1000




#%%
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Visualize fragment overlaps with peaks.")
    parser.add_argument("--res_dir", required=True, type=str, help="Directory to results.")
    
    args = parser.parse_args()
    output_dir = args.res_dir
    
    # Create subdirectory for fragments overlap peaks    
    figure_subdir = os.path.join(output_dir, 'figures', 'fragments_overlap_peaks')
    if not os.path.exists(figure_subdir):
        os.makedirs(figure_subdir)

    fragments_overlap_subdir = os.path.join(output_dir, 'fragments_overlap_peaks')


    total_frag_counts_file = os.path.join(fragments_overlap_subdir, 'total_fragments_counts.txt')
    overlap_peaks_counts_file = os.path.join(fragments_overlap_subdir, 'overlap_peaks_counts.txt')
    overlap_entropy_peaks_counts_file = os.path.join(fragments_overlap_subdir, 'overlap_entropy_peaks_counts.txt')

    total_frag_counts = pd.read_csv(total_frag_counts_file, sep=" ", index_col=1, header=None)
    overlap_peaks_counts = pd.read_csv(overlap_peaks_counts_file, sep=" ", index_col=1, header=None)
    overlap_entropy_peaks_counts = pd.read_csv(overlap_entropy_peaks_counts_file, sep=" ", index_col=1, header=None)


    total_frag_counts.columns = ['total_fragments']
    overlap_peaks_counts.columns = ['frag_overlap_peaks']
    overlap_entropy_peaks_counts.columns = ['frag_overlap_entropy_peaks']

    #%% 
    # Join fragments counts on barcode index
    frag_overlap_peaks_df = total_frag_counts.join(overlap_peaks_counts, how='outer')
    frag_overlap_entropypeaks_df = total_frag_counts.join(overlap_entropy_peaks_counts, how='outer')


    frag_overlap_peaks_df = frag_overlap_peaks_df.fillna(0)
    frag_overlap_entropypeaks_df = frag_overlap_entropypeaks_df.fillna(0)

    # calculate the percentage of fragments overlapping with peaks/entropy-peaks 
    frag_overlap_peaks_df['frag_overlap_peaks%'] = frag_overlap_peaks_df['frag_overlap_peaks'] / frag_overlap_peaks_df['total_fragments'] * 100
    frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'] = frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks'] / frag_overlap_entropypeaks_df['total_fragments'] * 100


    # Filter out barcodes with less than MininumFragments fragments
    overlap_peaks_df = frag_overlap_peaks_df[frag_overlap_peaks_df['total_fragments'] > MininumFragments].sort_values(by='frag_overlap_peaks%', ascending=False)
    overlap_entropy_peaks_df = frag_overlap_entropypeaks_df[frag_overlap_entropypeaks_df['total_fragments'] > MininumFragments].sort_values(by='frag_overlap_entropy_peaks%', ascending=False)

    #%%
    

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.sort(overlap_peaks_df['frag_overlap_peaks%'].values)[::-1], label='Overlap with Peaks', color='blue', marker='.', markersize=2, alpha=0.5)
    ax.plot(np.sort(overlap_entropy_peaks_df['frag_overlap_entropy_peaks%'].values)[::-1], label='Overlap with Entropy Peaks', color='orange', marker='.', markersize=2, alpha=0.5)
    ax.set_xlabel('Rank of barcodes ')
    ax.set_ylabel('Percentage of fragments')
    fig.legend()
    fig_file = os.path.join(figure_subdir, 'fragments_overlap_peaks.png')
    fig.savefig(fig_file, dpi=300, bbox_inches='tight')