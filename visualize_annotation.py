#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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


    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='red', alpha=0.5, label='Peaks')
    ax.plot(entropy_peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='blue', alpha=0.5, label='Entropy Peaks')
    ax.plot(lost_peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='gray', alpha=0.5, label='Lost Peaks')
    ax.plot(newly_discovered_peaks_annot['Peak Score'].sort_values(ascending=False), '-', color='green', alpha=0.5, label='Newly Discovered Peaks')
    ax.set_title('Peak Score Comparison')
    ax.set_xlabel('Peak rank')
    ax.set_ylabel('Peak Score')
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
    freq_df.plot(kind='bar', width=0.8, alpha=0.7, ax=ax)
    ax.set_title('Peak Annotation Frequency Comparison')
    ax.set_xlabel('Annotation')
    ax.set_ylabel('Frequency')
    ax.legend( )
    fig.tight_layout()
    fig_file = os.path.join(figures_subdir, 'peak_annotation_frequency.png')
    plt.savefig(fig_file)





