# Automatically find threshold for selecting barcodes based on entropy values

#%% 
# Set default parameters
_k = 1
_cutoff = 6e-8
_spline_s = 77
_spline_k = 3

#%%
import matplotlib.pyplot as plt
#import numpy as np
import pickle

import scipy
#from scipy.interpolate import BSpline, CubicSpline, make_interp_spline, make_splrep, splev, make_interp_spline
from typing import Iterable
import os


#%%
# Use the current directory for imports
from autothreshold_kevintopush import   fit_spline_and_find_cutoff,  get_x_y_from_pickle_helper,  make_plots


#%%
if __name__ == "__main__":
    # parerse args
    import argparse
    parser = argparse.ArgumentParser(description='Find threshold for selecting barcodes based on entropy values')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to output directory')
    parser.add_argument('--chromosome', type=str, required=True, default='entropies.pickle', help='Path to pickle file with entropy values')
    
    
    args = parser.parse_args()
    output_dir = args.output_dir
    chromosome = args.chromosome

    pickle_path = os.path.join(output_dir, f"{chromosome}_barcode_entropy.pickle")
    figure_subdir = os.path.join(output_dir, "figures")    
    
    # Entropy is at log10 scale
    rank_list, sorted_log10entropy_list = get_x_y_from_pickle_helper(pickle_path)
    
    # Find threshold using Spline fitting
    rank_cutoff, log10_entry_cutoff, cs = fit_spline_and_find_cutoff(rank_list, sorted_log10entropy_list, k=_k, cutoff=_cutoff, spline_s=_spline_s, spline_k=_spline_k)

    print(f"Rank cutoff: {rank_cutoff}, log10 Entropy cutoff: {log10_entry_cutoff}")

    with open(os.path.join(output_dir, f"entropy_cutoff.csv"), 'w') as f:
        f.write(f"entropy cutoff,{10**log10_entry_cutoff}\n")
        f.write(f"log10 entropy cutoff,{log10_entry_cutoff}\n")
        f.write(f"rank cutoff,{rank_cutoff}\n")
        
    
    # Make plots
    spl_y = cs(rank_list)
    deriv = {}
    for i in range(1, 3):
        deriv[i] = cs.derivative(i)(rank_list)

    fig, ax = plt.subplots(3, 1, figsize=(10,10), dpi=300)
    ax[0].plot(rank_list, sorted_log10entropy_list, 'o', label='data')
    ax[0].plot(rank_list, spl_y, label='spline')
    ax[1].plot(rank_list, deriv[1], label='first derivative')
    ax[2].plot(rank_list, deriv[2], label='second derivative')

    ax[1].set_ylim(-0.001, 0.001)
    ax[2].set_ylim(-1e-7, 1e-7)

    for i in range(3):
        ax[i].axvline(rank_cutoff, color='pink', linestyle='--', label='rank_cutoff')
        if i == 0:
            ax[i].axhline(log10_entry_cutoff, color='pink', linestyle='--', label='entry_cutoff')

    ax[0].set_xlabel('rank')
    ax[0].set_ylabel('log 10 entropy')
    ax[0].set_title('fit spline')
    ax[1].set_title('first derivative')
    ax[2].set_title('second derivative')

    fig.tight_layout()
    fig.savefig(os.path.join(figure_subdir, f"{chromosome}_entropy_threshold_fitting_k{_k}_s{_spline_s}.png"))
    plt.close(fig)