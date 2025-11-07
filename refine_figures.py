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


        total_frag_counts_file = os.path.join(output_dir, 'total_fragments_counts.txt')
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
    parser.add_argument('--entropythreshold', type=float, default=None, help='Entropy threshold for filtering.')
    parser.add_argument("--crbarcode_file", type=str, required=True, help="CellRanger barcode file, i.e., barcodes labeled to contain cells.")


    args = parser.parse_args()
    output_dir = args.output_dir
    chromosome = args.chromosome
    EntropyThreshold = args.entropythreshold
    crbarcode_file = args.crbarcode_file

    #%%
    import os 
    
    if EntropyThreshold is None:
        entropy_cutoff_file = os.path.join(output_dir, f"entropy_cutoff.csv")
        with open(entropy_cutoff_file, 'r') as f:
            line = f.readline()  # first line has entropy cutoff info
            if 'entropy cutoff' not in line:
                raise ValueError("File entropy_cutoff.csv format incorrect.")
            EntropyThreshold = float(line.strip('\n').split(',')[1])
    
    figure_subdir = os.path.join(output_dir, 'figures')

    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')
    
    with open(entropy_file, 'rb') as f:
        entropy_dict = pickle.load(f)
    
    # remove None values    
    noNone_entropies = {k: v for k, v in entropy_dict.items() if v is not None}

    
    entropy = np.array(list(noNone_entropies.values()))
    entropy.sort() # sort the entropy values


    # %%
    # Knee plot for entropy
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.log10(entropy[::-1] + precision), '-', color='blue', alpha=0.5, label='Entropy')
    ax.axhline(y=np.log10(EntropyThreshold + precision), color='red', linestyle='--', label=f'log10("{EntropyThreshold}")')
    ax.set_ylabel('Entropy (log10-scaled)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_logscale.png')
    fig.savefig(fig_file, bbox_inches='tight')
    #%%
    histogram of log10(entropy) values -- for (potentially) Gaussian mixture model fitting
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

    #%%
    Nbcs_to_plot = 30000  # 10x can only process up to 20k cells
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(bc_entropy_sort['index'][:Nbcs_to_plot], np.log10(bc_entropy_sort['entropy'][:Nbcs_to_plot] + precision), \
                c=bc_entropy_sort['cr_cell'][:Nbcs_to_plot].map({True: 'blue', False: 'orange'}), \
                marker= '.', alpha=0.4, label='CR Cell Barcodes', s=1)
    ax.axhline(y=np.log10(EntropyThreshold + precision), color='red', linestyle='--', label=f'log10({EntropyThreshold})')
    ax.set_ylabel('Entropy (log10-scaled)')
    ax.set_xlabel('Barcode Rank')
    fig.legend()
    fig_file = os.path.join(figure_subdir, f'{chromosome}_entropy_knee_plot_color_crbc_sub30k.png')
    fig.savefig(fig_file, bbox_inches='tight')
    

    
    
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
    
    #%%
    fig, ax = plt.subplots(figsize=(10, 6))
    p = ax.scatter(frag_overlap_entropypeaks_df['log10_entropy'], np.log10(frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks'] + precision), c=frag_overlap_entropypeaks_df['log10_total_fragments'], \
        cmap='viridis', alpha=0.6, s=10)
    fig.colorbar(p, ax=ax, label='log10(Total Fragments)')
    ax.set_xlabel('log10(Entropy)')
    ax.set_ylabel('log10(Fragment Overlap Entropy Peaks)')
    fig_file = os.path.join(figure_subdir, f'log10_entropy_vs_fragOverlapEntropyPeaks.png')
    fig.savefig(fig_file)
    

    #%% 
    # Color by CellRanger cell-calling
    frag_overlap_entropypeaks_df['cr_cell'] = False
    frag_overlap_entropypeaks_df.loc[frag_overlap_entropypeaks_df.index.isin(cr_bc[0]), 'cr_cell'] = True
    fig, ax = plt.subplots(figsize=(10, 6))
    p = ax.scatter(frag_overlap_entropypeaks_df[ 'log10_entropy'], 
                   np.log10(frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks'] + precision), \
                   c=frag_overlap_entropypeaks_df['cr_cell'], alpha=0.6, s=10)
    ax.set_xlabel('log10(Entropy)')
    ax.set_ylabel('log10(Fragment Overlap Entropy Peaks)')
    ax.legend(*p.legend_elements(), title="CR Cell", loc='lower right')
    fig_file = os.path.join(figure_subdir, f'log10_entropy_vs_fragOverlapEntropyPeaks_CRCalling.png')
    fig.savefig(fig_file)


    fig, ax = plt.subplots(figsize=(10, 6))
    p = ax.scatter(frag_overlap_entropypeaks_df['log10_total_fragments'], frag_overlap_entropypeaks_df['log10_entropy'], c=frag_overlap_entropypeaks_df['cr_cell'], \
        cmap='viridis', alpha=0.6, s=10)
    ax.legend(*p.legend_elements(), title="CR Cell", loc='lower right')
    ax.set_ylabel('log10(Entropy)')
    ax.set_xlabel('log10(Total Fragments)')
    fig_file = os.path.join(figure_subdir, f'log10_entropy_vs_total_fragments_CRCalling.png')
    fig.savefig(fig_file)


# %%
#### TODO 
    frag_df = pd.read_csv(frag_file, sep='\t', header=None)
    frag_df.columns = ['chrom', 'left_insert', 'right_insert', 'cb', 'support']
    frag_df['length'] = frag_df['right_insert'] - frag_df['left_insert']
    filter_frag_df = frag_df[frag_df['cb'].isin(bc_pass_entropy.index)]
    garbage_frag_df = frag_df[~ frag_df['cb'].isin(filter_frag_df['cb'])]

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), sharey=True)
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
    fig_file = os.path.join(output_dir, 'figures', 'frag_df_length.png')
    fig.tight_layout()
    fig.savefig(fig_file)
