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
precision = 1e-10
# if __name__ == "__main__":
# import argparse
# parser = argparse.ArgumentParser(description="Overlay RNA cell type annotation on ATAC peak latent space")
# parser.add_argument('--res_dir', type=str, required=True, help='Directory containing the results')
# args = parser.parse_args()

# output_dir = args.res_dir
output_dir = '/home/syyang/adipose_ln/atac/res/brain2_3k'
chromosome = 'chr1'

#%%
# load fragments overlap peaks data & entropy data 
frag_overlap_peaks_df, frag_overlap_entropypeaks_df = get_frag_overlap_peaks_df(output_dir)


DownstreamReanalysis_dir = os.path.join(output_dir, 'DownstreamReanalysis')
bc_peak_count_h5 = os.path.join(DownstreamReanalysis_dir, 'peakvi.h5ad')
entropy_bc_peak_count_h5 = os.path.join(DownstreamReanalysis_dir, 'entropy_peakvi.h5ad')

#%%
peakvi_adata = sc.read_h5ad(bc_peak_count_h5)
peakvi_adata.obs['atac_barcodes'] = peakvi_adata.obs['barcodes'].apply(lambda x: x.split('-')[0])

entropy_peakvi_adata = sc.read_h5ad(entropy_bc_peak_count_h5)
entropy_peakvi_adata.obs['atac_barcodes'] = entropy_peakvi_adata.obs['barcodes'].apply(lambda x: x.split('-')[0])

# %%
#%%
# Since brain2_3k was processed using cellranger-arc. (tho I only used the atac-fragments)
# the barcodes was already mapped to the rna barcodes.
# So there is no need to map again.
# Hence commenting out the following lines:
# ## load rna and atac barcode map 
# atac_rna_barcode_file = '/mnt/files/syy/adipose/multiom/atac_rna_barcodes_map.tsv'
# atac_rna_barcode_map = pd.read_csv(atac_rna_barcode_file, sep='\t')
# peakvi_adata.obs = peakvi_adata.obs.merge(atac_rna_barcode_map, left_on='atac_barcodes', right_on='atac_barcodes')
# entropy_peakvi_adata.obs = entropy_peakvi_adata.obs.merge(atac_rna_barcode_map, left_on='atac_barcodes', right_on='atac_barcodes')

#%%
## load cell type annotation analyzed using RNA only 
rna_celltype_annotation_file = '/mnt/hdd_bob/syy/adipose/atac/cr_goldstand/brain2_3k/analysis/clustering/gex/graphclust/clusters.csv'
rna_celltype_annotation = pd.read_csv(rna_celltype_annotation_file,  header=0)
rna_celltype_annotation['rna_cluster'] = rna_celltype_annotation['Cluster'].astype(str)
# %%


#%%
peakvi_adata.obs = peakvi_adata.obs.merge(rna_celltype_annotation[['Barcode', 'rna_cluster']],
                    left_on='barcodes', right_on='Barcode', how='left') 
entropy_peakvi_adata.obs = entropy_peakvi_adata.obs.merge(rna_celltype_annotation[['Barcode', 'rna_cluster']],
                    left_on='barcodes', right_on='Barcode', how='left')
# %%

figure_dir = os.path.join(DownstreamReanalysis_dir, 'figures')
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)


# save the adata with cell type annotation
os.chdir(DownstreamReanalysis_dir)
peakvi_fig_file = os.path.join(figure_dir, 'umap_peakvi.png')
sc.pl.umap(peakvi_adata, color=['rna_cluster', 'clusters_peakvi'], ncols=1,  title='RNA-cell-type on atac peaks latent space', \
        size=30, na_color='None', save='umap_peakvi.png')
#%%
entropy_peakvi_fig_file = os.path.join(figure_dir, 'umap_entropy_peakvi.png')
sc.pl.umap(entropy_peakvi_adata, color=['rna_cluster', 'clusters_peakvi'], ncols=1,  title='RNA-cell-type on atac entropy peaks latent space',  \
        size=30, na_color='None', save='umap_entropy_peakvi.png')

# %%

peakvi_adata.obs = peakvi_adata.obs.set_index('barcodes').join(frag_overlap_peaks_df, how='left')
peakvi_adata.obs['log10_entropy'] = np.log10(peakvi_adata.obs['entropy'] + precision)

#%%
entropy_peakvi_adata.obs = entropy_peakvi_adata.obs.set_index('barcodes').join(frag_overlap_entropypeaks_df, how='left')
entropy_peakvi_adata.obs['log10_entropy'] = np.log10(entropy_peakvi_adata.obs['entropy'] + precision)

#%%
peakvi_adata.obs['log2_total_fragments'] = np.log2(peakvi_adata.obs['total_fragments'] + precision)
entropy_peakvi_adata.obs['log2_total_fragments'] = np.log2(entropy_peakvi_adata.obs['total_fragments'] + precision)

#%%
sc.pl.umap(peakvi_adata, color=['frag_overlap_peaks%', 'log10_entropy', 'log2_total_fragments'], ncols=3, 
        title='Percentage of fragments overlapping with peaks and entropy', size=30, na_color='None', 
        save='umap_peakvi_bc_qc.png')

#%%
sc.pl.umap(entropy_peakvi_adata, color=['frag_overlap_entropy_peaks%', 'log10_entropy', 'log2_total_fragments'], ncols=3,
        title='Percentage of fragments overlapping with entropy peaks and entropy', size=30, na_color='None', 
        save='umap_entropy_peakvi_bc_qc.png')


# %%
