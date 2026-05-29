#%%
import matplotlib.pyplot as plt
import numpy as np

#%%
import sys
import os

# Get the absolute path to folder-A (the parent of the current script's folder)
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from autothreshold import   fit_spline_and_find_cutoff,  get_sorted_entropy_helper

#%%
windowsize_list='100 500 1000 3000 6000 10000 15000 20000 30000 50000 80000 100000 200000 400000 800000'
ws_list = list(map(int, windowsize_list.split()))

entropy_ws_pdir = '/mnt/hdd_bob/syy/adipose/atac/res'

sample_list = ['VIB_10xmultiome_2', 'mouse_cryo_brain_5k', 'human_10k_1', 'mice_brain_5k']
Ncells_list = {}

# sample='VIB_10xmultiome_2'
for sample in sample_list:
    
    Ncells_list[sample] = []
    for ws in ws_list:
        entropy_dir = os.path.join(entropy_ws_pdir, f'_{sample}_WS{ws}')
        entropy_cutoff_file = os.path.join(entropy_dir, 'entropy_cutoff.csv')
        if not os.path.exists(entropy_cutoff_file):
            print(f'Entropy cutoff file not found for window size {ws}. Skipping.')
            continue
        with open(entropy_cutoff_file, 'r') as f:
            
            cutoff_value = f.read().split('\n')[-2].split(',')
            print(f'{sample}, Window size: {ws}, {cutoff_value}')
            Ncells_list[sample].append(int(int(cutoff_value[1])))

#%%
_k = 3
_spline_s = 1
_spline_k = 3

output_dir = os.path.join(entropy_ws_pdir, f'_human_10k_1_WS6000')
pickle_path = os.path.join(output_dir, f"calculated_barcode_entropy.pickle")
sorted_log10entropy_list = get_sorted_entropy_helper(pickle_path)

rank_cutoff, log10_entry_cutoff, cs, inflection_regions, regions_rank = fit_spline_and_find_cutoff(sorted_log10entropy_list, k=_k, spline_s=_spline_s, spline_k=_spline_k)

print(rank_cutoff)
# %%


# Make plots
rank_list = np.arange(len(sorted_log10entropy_list))
spl_y = cs(rank_list)
deriv = {}
for i in range(1, 3):
    deriv[i] = cs.derivative(i)(rank_list)

deriv[2] = deriv[2] / np.abs(deriv[2]).max()  # normalize second derivative for better visualization    

fig, ax = plt.subplots(3, 1, figsize=(8,6), dpi=800)
ax[0].plot(rank_list, sorted_log10entropy_list, 'o', label='data')
ax[0].plot(rank_list, spl_y, label='spline')
ax[1].plot(rank_list, deriv[1], label='first derivative')
ax[2].plot(rank_list, deriv[2], label='second derivative')

d1_limit = np.abs(deriv[1]).max()
ax[1].set_ylim(-d1_limit*1.1, d1_limit*1.1)
# %%
