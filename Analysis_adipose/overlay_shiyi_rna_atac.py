#%%
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


        frag_overlap_peaks_df = frag_overlap_peaks_df.fillna(0)
        frag_overlap_entropypeaks_df = frag_overlap_entropypeaks_df.fillna(0)

        # calculate the percentage of fragments overlapping with peaks/entropy-peaks 
        frag_overlap_peaks_df['frag_overlap_peaks%'] = frag_overlap_peaks_df['frag_overlap_peaks'] / frag_overlap_peaks_df['total_fragments'] * 100
        frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks%'] = frag_overlap_entropypeaks_df['frag_overlap_entropy_peaks'] / frag_overlap_entropypeaks_df['total_fragments'] * 100

        # Join entropy data on barcode index
        frag_overlap_peaks_df = frag_overlap_peaks_df.join(entropy_df, how='left')
        frag_overlap_entropypeaks_df = frag_overlap_entropypeaks_df.join(entropy_df, how='left')

        return [frag_overlap_peaks_df, frag_overlap_entropypeaks_df]
#%%

import scanpy as sc
import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt

#%%
MininumFragments = 1000
precision = 1e-10
# if __name__ == "__main__":
# import argparse
# parser = argparse.ArgumentParser(description="Overlay RNA cell type annotation on ATAC peak latent space")
# parser.add_argument('--res_dir', type=str, required=True, help='Directory containing the results')
# parser.add_argument('--rna_analysis_dir', type=str, required=True, help='Directory containing the RNA results')
# args = parser.parse_args()

# output_dir = args.res_dir
# corresponding_rna_dir = args.rna_analysis_dir
output_dir = '/home/syyang/adipose_ln/atac/res/NK_ATAC_MAPLE013_CS'
corresponding_rna_dir = '/mnt/hdd_alice/syy/adipose/rna/Analysis/QC_Process_output/NK_GEX_MAPLE013_CS'

chromosome = 'chr1'


#%%
# load fragments overlap peaks data & entropy data 
frag_overlap_peaks_df, frag_overlap_entropypeaks_df = get_frag_overlap_peaks_df(output_dir)


#%%
sample_name = os.path.basename(output_dir)
rna_sample_name = sample_name.split('ATAC_')[1]

DownstreamReanalysis_dir = os.path.join(output_dir, 'DownstreamReanalysis')
bc_peak_count_h5 = os.path.join(DownstreamReanalysis_dir, 'peakvi.h5ad')
entropy_bc_peak_count_h5 = os.path.join(DownstreamReanalysis_dir, 'entropy_peakvi.h5ad')

#%%
peakvi_adata = sc.read_h5ad(bc_peak_count_h5)
peakvi_adata.obs['atac_barcodes'] = peakvi_adata.obs['barcodes'].apply(lambda x: x.split('-')[0])

entropy_peakvi_adata = sc.read_h5ad(entropy_bc_peak_count_h5)
entropy_peakvi_adata.obs['atac_barcodes'] = entropy_peakvi_adata.obs['barcodes'].apply(lambda x: x.split('-')[0])

# %%
## load rna and atac barcode map 
atac_rna_barcode_file = '/mnt/files/syy/adipose/multiom/atac_rna_barcodes_map.tsv'
atac_rna_barcode_map = pd.read_csv(atac_rna_barcode_file, sep='\t')

#%%

peakvi_adata.obs = peakvi_adata.obs.merge(atac_rna_barcode_map, left_on='atac_barcodes', right_on='atac_barcodes')
entropy_peakvi_adata.obs = entropy_peakvi_adata.obs.merge(atac_rna_barcode_map, left_on='atac_barcodes', right_on='atac_barcodes')

#%%
## load cell type annotation analyzed using RNA only 
rna_celltype_annotation_file = '/home/syyang/GitRepo/atac/Analysis_adipose/allsamples_rna/allsample_rna_celltype_asign_sean.tsv'
rna_celltype_annotation = pd.read_csv(rna_celltype_annotation_file, sep='\t', index_col=0)
rna_celltype_annotation['Sample'] = rna_celltype_annotation['sample'].apply(lambda x: x.upper()) # make it sample name upper case
rna_celltype_annotation['rna_barcodes'] = rna_celltype_annotation.apply(lambda x: x.name.split('-')[0], axis = 1)
# %%

sample_rna_celltype = rna_celltype_annotation[rna_celltype_annotation['Sample'] == sample_name.split('ATAC_')[0] + 'GEX_' +rna_sample_name].copy()

#%%
peakvi_adata.obs = peakvi_adata.obs.merge(sample_rna_celltype[['rna_barcodes', 'cell_type']],
                    left_on='rna_barcodes', right_on='rna_barcodes', how='left') 
entropy_peakvi_adata.obs = entropy_peakvi_adata.obs.merge(sample_rna_celltype[['rna_barcodes', 'cell_type']],
                    left_on='rna_barcodes', right_on='rna_barcodes', how='left')
# %%

figure_dir = os.path.join(DownstreamReanalysis_dir, 'figures')
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)


# save the adata with cell type annotation
os.chdir(DownstreamReanalysis_dir)
sc.pl.umap(peakvi_adata, color=['cell_type', 'clusters_peakvi'], ncols=1,  title='RNA-cell-type on atac peaks latent space', \
        size=30, na_color='None'
        , save='umap_peakvi.png')
sc.pl.umap(entropy_peakvi_adata, color=['cell_type', 'clusters_peakvi'], ncols=1,  title='RNA-cell-type on atac entropy peaks latent space',  \
        size=30, na_color='None'
        , save='umap_entropy_peakvi.png')

# %%
# Overlay with entropy and fragments overlap peaks percentage



# %%
peakvi_adata.obs = peakvi_adata.obs.set_index('barcodes').join(frag_overlap_peaks_df, how='left')
peakvi_adata.obs['log10_entropy'] = np.log10(peakvi_adata.obs['entropy'] + precision)

#%%
entropy_peakvi_adata.obs = entropy_peakvi_adata.obs.set_index('barcodes').join(frag_overlap_entropypeaks_df, how='left')
entropy_peakvi_adata.obs['log10_entropy'] = np.log10(entropy_peakvi_adata.obs['entropy'] + precision)

# %%
# Project fragments% overlap peaks and entropy on umap
sc.pl.umap(peakvi_adata, color=['frag_overlap_peaks%', 'log10_entropy'], ncols=2, 
        title='Percentage of fragments overlapping with peaks and entropy', size=30, na_color='None', 
        save='umap_peakvi_bc_qc.png')

sc.pl.umap(entropy_peakvi_adata, color=['frag_overlap_entropy_peaks%', 'log10_entropy'], ncols=2,
        title='Percentage of fragments overlapping with entropy peaks and entropy', size=30, na_color='None', 
        save='umap_entropy_peakvi_bc_qc.png')

# %%


# Load RNA doublet scores
rna_doublet_file = os.path.join(corresponding_rna_dir, 'scds_doublets_singlets.tsv')
rna_doublet_df = pd.read_csv(rna_doublet_file, sep='\t', index_col=0)
rna_doublet_df['rna_barcodes'] = rna_doublet_df.apply(lambda x: x.name.split('-')[0], axis = 1)

#
peakvi_adata.obs = peakvi_adata.obs.merge(rna_doublet_df[['rna_barcodes', 'scds_score', 'scds_DropletType']], left_on='rna_barcodes', right_on='rna_barcodes', how='left')
entropy_peakvi_adata.obs = entropy_peakvi_adata.obs.merge(rna_doublet_df[['rna_barcodes', 'scds_score', 'scds_DropletType']], left_on='rna_barcodes', right_on='rna_barcodes', how='left')
#%%
# Project doublet scores on umap

sc.pl.umap(entropy_peakvi_adata, color=['scds_score', 'scds_DropletType', 'total_fragments', 'log10_entropy'], ncols=2,
        title='doublet score on entropy peak atac', size=30, na_color='None',
        save='umap_entropy_peakvi_overlay_rnadoublet.png')
# %%
fig_file = os.path.join(figure_dir, 'violin_entropy_rnadoublet.png')
fig, ax = plt.subplots(ncols=   2, sharey=True)
ax[0].violinplot(peakvi_adata.obs[peakvi_adata.obs['scds_DropletType'] == 'doublet']['log10_entropy'], showmeans=True, showmedians=True, bw_method=0.1)
ax[1].violinplot(peakvi_adata.obs[~ (peakvi_adata.obs['scds_DropletType'] == 'doublet')]['log10_entropy'], showmeans=True, showmedians=True, bw_method=0.1)
ax[0].set_xlabel('doublet')
ax[1].set_xlabel('singlet')
ax[0].set_ylabel('log10 Entropy')
fig.legend()
fig.savefig(fig_file)
# %%
