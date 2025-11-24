#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import get_frag_overlap_peaks_df


precision = 1e-10

#%%
def simplify_annotation(annot_df, inplace=False):
    """
    Simplify the Annotation column.
    """
    if not inplace:
        annot_df = annot_df.copy()
        annot_df['simple_annotation'] = annot_df['Annotation'].apply(lambda x: str(x).split('(')[0])  
        return annot_df
    else:
        # in place change of the  DataFrame
        annot_df['simple_annotation'] = annot_df['Annotation'].apply(lambda x: str(x).split('(')[0])  
        



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Visualize ATAC-seq annotation')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for the results')

    args = parser.parse_args()
    output_dir = args.output_dir

    #%% 
    import os
    annotation_dir = os.path.join(output_dir, 'annotation')

    figures_subdir = os.path.join(output_dir, 'figures', 'peak_annotation')
    if not os.path.exists(figures_subdir):
        os.makedirs(figures_subdir)
    # %%

    peaks_annot = pd.read_csv(os.path.join(annotation_dir, 'annotated_peaks.bed'), \
        sep='\t')
    # %%
    entropy_peaks_annot = pd.read_csv(os.path.join(annotation_dir, 'annotated_entropy_peaks.bed'), \
        sep='\t')
    # %%

    lost_peaks_annot = pd.read_csv(os.path.join(annotation_dir, 'annotated_lost_peaks.bed'), \
        sep='\t')
    # %%
    newly_discovered_peaks_annot = pd.read_csv(os.path.join(annotation_dir, 'annotated_newly_discovered_peaks.bed'), \
        sep='\t')
    #%% 

    simplify_annotation(peaks_annot, inplace=True)
    simplify_annotation(entropy_peaks_annot, inplace=True)
    # %%


    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='orange', alpha=0.5, label='Peaks', lw=4)
    ax.plot(entropy_peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='blue', alpha=0.5, label='Entropy Peaks', lw=4)
    ax.plot(lost_peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='gray', alpha=0.5, label='Lost Peaks', lw=4)
    ax.plot(newly_discovered_peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='green', alpha=0.5, label='Newly Discovered Peaks', lw=4)
    ax.set_title('Peak Score Comparison', fontsize=14)
    ax.set_xlabel('Peak rank', fontsize=14)
    ax.set_ylabel('Peak Score', fontsize=14)
    ax.legend( loc='upper right', bbox_to_anchor=(1.6, .8))
    fig.tight_layout()
    fig_file = os.path.join(figures_subdir, 'peak_score_comparison.png')
    plt.savefig(fig_file)
    # %%



    #%% 
    peaks_freq = peaks_annot['simple_annotation'].value_counts()
    entropy_peaks_freq = entropy_peaks_annot['simple_annotation'].value_counts()
    freq_df = pd.DataFrame({
        'Peaks': peaks_freq,
        'Entropy Peaks': entropy_peaks_freq
    }).fillna(0)

    fig, ax = plt.subplots(figsize=(10, 6))
    freq_df.plot(kind='bar', width=0.8, alpha=0.9, ax=ax, color=['orange', 'blue'], fontsize=14)
    ax.set_title('Peak Annotation Frequency Comparison', fontsize=16)
    ax.set_xlabel('Annotation', fontsize=16)
    ax.set_ylabel('Frequency', fontsize=16)
    ax.legend()
    fig.tight_layout()
    fig_file = os.path.join(figures_subdir, 'peak_annotation_frequency.png')
    plt.savefig(fig_file)


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
    p = ax.scatter(frag_overlap_entropypeaks_df['log10_entropy'], np.log10(frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks'] + precision), c=frag_overlap_entropypeaks_df['log10_total_fragments'], \
        cmap='viridis', alpha=0.6, s=10)
    fig.colorbar(p, ax=ax, label='log10(Total Fragments)')
    ax.set_xlabel('log10(Entropy)')
    ax.set_ylabel('log10(Fragment Overlap Entropy Peaks)')
    fig_file = os.path.join(figure_subdir, f'log10_entropy_vs_fragOverlapEntropyPeaks.png')
    fig.savefig(fig_file)






