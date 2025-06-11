import matplotlib.pyplot as plt
import pickle 

#### GLOBAL VARIABLES >>
EntropyThreshold = 0.000609


#### GLOBAL VARIABLES << 


# %%
#  >>>>>>>>>..
# >>>>>>>>>>..
fl_entropy = pd.DataFrame.from_dict(fl_entropy, orient='index', columns=['entropy'])
# %%

fig, ax = plt.subplots()
ax.hist(fl_entropy['entropy'], bins=100, alpha=0.5, label='fl', color='green')
fig.legend()
# %%

## see if barcodes filtering in CR would separate the entropy distributions within each sample
fl_bc = pd.read_csv('/mnt/files/syy/adipose/atac/cr_atac_results/AB_ATAC_FL_TAM_SQ/outs/filtered_peak_bc_matrix/barcodes.tsv', sep='\t', header=None)

#%%


fl_bc_entropy = fl_entropy.merge(fl_bc, left_index=True, right_on=0, how='right')
fl_bc_entropy = fl_bc_entropy.set_index(0)


fl_empty_entropy = fl_entropy[fl_entropy.index.isin(fl_bc[0]) == False]


# %%
# %%
fig, ax = plt.subplots(1, 2, sharey=True)
ax[0].violinplot(fl_bc_entropy['entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
ax[1].violinplot(fl_empty_entropy['entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
ax[0].set_ylabel('Entropy')
ax[0].set_xlabel('CR cell barcodes')
ax[1].set_xlabel('CR empty barcodes')
#0.000609 is 75% from the fl empty barcodes
ax[0].axhline(y=EntropyThreshold, color='r', linestyle='-', label=f'{EntropyThreshold}')
ax[1].axhline(y=EntropyThreshold, color='r', linestyle='-', label=f'{EntropyThreshold}')
fig.legend()
fig.suptitle('AB_ATAC_FL_TAM_SQ dataset')
fig.savefig('fl_entropy.png')
# %%

# now provide only BC that pass entropy threshold
fl_entropy['pass'] = fl_entropy['entropy'] > EntropyThreshold



# %%
fl_entropy.loc[fl_entropy['pass'] == True].to_csv(os.path.join(res_dir, 'fl_pass_entropy.tsv'), sep='\t')
# %%

# load filtered fragments based on bc thresholded above 
filtered_frag_file = 'pass_entropy_bc/fl_fragments.tsv'
fl_frag_file= '/mnt/files/syy/adipose/atac/cr_atac_results/eda_AB_ATAC_FL_TAM_SQ/fragments.tsv'

fl_frag = pd.read_csv(frag, sep='\t', header=None, comment='#')
fl_filter_frag = pd.read_csv(os.path.join(res_dir, filtered_frag_file), sep='\t', header=None)

#%%
fl_frag.columns = ['chrom', 'left_insert', 'right_insert', 'cb', 'support']
fl_frag['length'] = fl_frag['right_insert'] - fl_frag['left_insert']
fl_filter_frag.columns = ['chrom', 'left_insert', 'right_insert', 'cb', 'support']
fl_filter_frag['length'] = fl_filter_frag['right_insert'] - fl_filter_frag['left_insert']
# %%

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6), sharey=True)
ax[0].hist(fl_filter_frag['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
ax[0].set_xlabel('Fragment length')
ax[0].set_ylabel('Frequency')
ax[0].set_title('Fragment length distribution for filtered fragments')
ax[1].hist(fl_frag['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
ax[1].set_xlabel('Fragment length')
ax[1].set_ylabel('Frequency')
ax[1].set_title('Fragment length distribution for all fragments')
fig.savefig('fl_frag_length.png')

