## Universal functions


import os
import numpy as np
import pandas as pd



# For defaultdict default factory 
# necessary for loading pickle files of Tn5_insertion_count results
def return_none():
    return None




def get_frag_overlap_peaks_df(output_dir, chromosome='chr1'):
        import pickle
        # Load entropy data
        entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')

        # subdirectory for fragments overlap peaks    
        fragments_overlap_subdir = os.path.join(output_dir, 'fragments_overlap_peaks')

        total_frag_counts_file = os.path.join(output_dir, 'total_fragments_counts.txt')
        overlap_peaks_counts_file = os.path.join(fragments_overlap_subdir, 'overlap_peaks_counts.txt')
        overlap_entropy_peaks_counts_file = os.path.join(fragments_overlap_subdir, 'overlap_entropy_peaks_counts.txt')


        # Validate all files exist
        required_files = {
                'entropy_file': entropy_file,
                'total_frag_counts': total_frag_counts_file,
                'overlap_peaks_counts': overlap_peaks_counts_file,
                'overlap_entropy_peaks_counts': overlap_entropy_peaks_counts_file
        }
        for name, filepath in required_files.items():
                if not os.path.exists(filepath):
                        raise FileNotFoundError(
                                f"Required file not found: {filepath}\n"
                                f"Please ensure the previous steps completed successfully."
                        )

        # Load entropy data
        with open(entropy_file, 'rb') as f:
                entropy_data = pickle.load(f)
        entropy_df = pd.DataFrame.from_dict(entropy_data, orient='index', columns=['entropy'])
        
        # Load fragment counts
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