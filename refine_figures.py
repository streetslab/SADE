#%%
import pandas as pd


def get_frag_overlap_peaks_df(output_dir, chromosome='chr1'):
        import pickle
        # Load entropy data
        entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')

        with open(entropy_file, 'rb') as f:
                entropy_data = pickle.load(f)

        
        entropy_df = pd.DataFrame.from_dict(entropy_data, orient='index', columns=['entropy'])

         
        # subdirectory for fragments overlap peaks    
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
        
        # Join fragments counts on barcode index
        frag_overlap_peaks_df = total_frag_counts.join(overlap_peaks_counts, how='outer')
        frag_overlap_entropypeaks_df = total_frag_counts.join(overlap_entropy_peaks_counts, how='outer')

        frag_overlap_peaks_df = frag_overlap_peaks_df.mask(frag_overlap_peaks_df.isna(), 0.)
        frag_overlap_entropypeaks_df = frag_overlap_entropypeaks_df.mask(frag_overlap_entropypeaks_df.isna(), 0.)


        # calculate the percentage of fragments overlapping with peaks/entropy-peaks 
        frag_overlap_peaks_df['frag_overlap_peaks%'] = frag_overlap_peaks_df['frag_overlap_peaks'] / frag_overlap_peaks_df['total_fragments'] * 100
        frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'] = frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks'] / frag_overlap_entropypeaks_df['total_fragments'] * 100

        # Join entropy data on barcode index
        frag_overlap_peaks_df = frag_overlap_peaks_df.join(entropy_df, how='left')
        frag_overlap_entropypeaks_df = frag_overlap_entropypeaks_df.join(entropy_df, how='left')

        return [frag_overlap_peaks_df, frag_overlap_entropypeaks_df]
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
    ### >>>> This chunk of code could be deleted (in calculate_entropy.py file too)
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
    ### <<<< This chunk of code could be deleted (in calculate_entropy.py file too)

    # %%
    
    bc_entropy = pd.DataFrame.from_dict(entropy_dict, orient='index', columns=['entropy'])
    
    cr_bc = pd.read_csv(crbarcode_file, sep='\t', header=None)
    cr_bc_pd = cr_bc.set_index(0)
    cr_bc_pd['cr_cell'] = True

    bc_entropy_sort = bc_entropy.sort_values(by='entropy', ascending=False)
    bc_entropy_sort['index'] = np.arange(bc_entropy_sort.shape[0])
    bc_entropy_sort = bc_entropy_sort.join(cr_bc_pd, how='left')
    bc_entropy_sort['cr_cell'] = bc_entropy_sort['cr_cell'].mask(bc_entropy_sort['cr_cell'].isna(), False)

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
    
    
    ### >>>> This chunk of code could be deleted (in overlay_entropy_CRcelllabel_plot.py file too)
    # Violin plot of entropy values for CellRanger labeled cells vs empty barcodes
    cr_bc_entropy = bc_entropy.merge(cr_bc, left_index=True, right_on=0, how='right')
    cr_bc_entropy = cr_bc_entropy.set_index(0)

    cr_empty_entropy = bc_entropy[bc_entropy.index.isin(cr_bc[0]) == False].copy()
    
    cr_bc_entropy['log10_entropy'] = np.log10(cr_bc_entropy['entropy'] + precision)
    cr_empty_entropy['log10_entropy'] = np.log10(cr_empty_entropy['entropy'] + precision)
    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].violinplot(cr_bc_entropy['log10_entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
    ax[1].violinplot(cr_empty_entropy['log10_entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
    ax[0].set_ylabel('log10 Entropy')
    ax[0].set_xlabel('CR cell barcodes')
    ax[1].set_xlabel('CR empty barcodes')
    ax[0].axhline(y=np.log10(EntropyThreshold + precision), color='r', linestyle='-', label=f'log10({EntropyThreshold})')
    ax[1].axhline(y=np.log10(EntropyThreshold + precision), color='r', linestyle='-', label=f'log10({EntropyThreshold})')
    fig.legend()
    fig_file = os.path.join(figure_subdir, 'log10_entropy_violin_EmptyvsCell.png')
    fig.savefig(fig_file)
    ### <<<< This chunk of code could be deleted (in overlay_entropy_CRcelllabel_plot.py file too)
    
    
    #%%
    # Scatter plot of log10(entropy) vs log10(total fragments) colored by frag% overlap peaks
    # load fragments overlap peaks data & entropy data 
    _, frag_overlap_entropypeaks_df = get_frag_overlap_peaks_df(output_dir)

    frag_overlap_entropypeaks_df['log10_total_fragments'] = np.log10(frag_overlap_entropypeaks_df['total_fragments'] + precision)
    frag_overlap_entropypeaks_df['log10_entropy'] = np.log10(frag_overlap_entropypeaks_df['entropy'] + precision)
    fig, ax = plt.subplots(figsize=(10, 6))
    p = ax.scatter(frag_overlap_entropypeaks_df['log10_total_fragments'], frag_overlap_entropypeaks_df['log10_entropy'], c=frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'], \
        cmap='viridis', alpha=0.6, s=10)
    fig.colorbar(p, ax=ax, label='Fragment Overlap Entropy Peaks %')
    ax.set_ylabel('log10(Entropy)')
    ax.set_xlabel('log10(Total Fragments)')
    fig_file = os.path.join(figure_subdir, f'log10_entropy_vs_total_fragments.png')
    fig.savefig(fig_file)
    
    #%%
    fig, ax = plt.subplots(figsize=(10, 6))
    p = ax.scatter(frag_overlap_entropypeaks_df['entropy'], frag_overlap_entropypeaks_df['total_fragments'], c=frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'], \
        cmap='viridis', alpha=0.6, s=10)
    fig.colorbar(p, ax=ax, label='Fragment Overlap Entropy Peaks %')
    ax.set_xlabel('Entropy')
    ax.set_ylabel('Total Fragments')
    fig_file = os.path.join(figure_subdir, f'entropy_vs_total_fragments.png')
    fig.savefig(fig_file)

# %%
