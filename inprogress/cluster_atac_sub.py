#%%

import scanpy as sc
import numpy as np
import pandas as pd
import os 
import scvi 
#%%


output_dir = '/mnt/hdd_bob/syy/adipose/atac/res/NK_ATAC_MAPLE013_SC'
sample_name = os.path.basename(output_dir)
rna_sample_name = sample_name.split('ATAC_')[1]

DownstreamReanalysis_dir = os.path.join(output_dir, 'DownstreamReanalysis')
bc_peak_count_h5 = os.path.join(DownstreamReanalysis_dir, 'bc_peak_count.h5ad')
entropy_bc_peak_count_h5 = os.path.join(DownstreamReanalysis_dir, 'entropy_bc_peak_count.h5ad')

#%%
bc_peak_count_adata = sc.read_h5ad(bc_peak_count_h5)
bc_peak_count_adata.obs['atac_barcodes'] = bc_peak_count_adata.obs['barcodes'].apply(lambda x: x.split('-')[0])

entropy_bc_peak_count_adata = sc.read_h5ad(entropy_bc_peak_count_h5)
entropy_bc_peak_count_adata.obs['atac_barcodes'] = entropy_bc_peak_count_adata.obs['barcodes'].apply(lambda x: x.split('-')[0])

# %%
## load rna and atac barcode map 
atac_rna_barcode_file = '/mnt/files/syy/adipose/multiom/atac_rna_barcodes_map.tsv'
atac_rna_barcode_map = pd.read_csv(atac_rna_barcode_file, sep='\t')

#%%

bc_peak_count_adata.obs = bc_peak_count_adata.obs.merge(atac_rna_barcode_map, left_on='atac_barcodes', right_on='atac_barcodes')
entropy_bc_peak_count_adata.obs = entropy_bc_peak_count_adata.obs.merge(atac_rna_barcode_map, left_on='atac_barcodes', right_on='atac_barcodes')

#%%
## load cell type annotation analyzed using RNA only 
rna_celltype_annotation_file = '/home/syyang/GitRepo/atac/Analysis_adipose/Henrique_Sean_RNAonly_celltype_multiome_samples.csv'
rna_celltype_annotation = pd.read_csv(rna_celltype_annotation_file, sep=',', index_col=0)
rna_celltype_annotation['Sample'] = rna_celltype_annotation['sample'].apply(lambda x: x.upper()) # make it sample name upper case
rna_celltype_annotation['rna_barcodes'] = rna_celltype_annotation['barcode'].apply(lambda x: x.split('-')[0])
# %%

sample_rna_celltype = rna_celltype_annotation[rna_celltype_annotation['Sample'] == rna_sample_name].copy()

#%%
bc_peak_count_adata.obs = bc_peak_count_adata.obs.merge(sample_rna_celltype[['rna_barcodes', 'cell_type']],
                       left_on='rna_barcodes', right_on='rna_barcodes', how='left') 
entropy_bc_peak_count_adata.obs = entropy_bc_peak_count_adata.obs.merge(sample_rna_celltype[['rna_barcodes', 'cell_type']],
                       left_on='rna_barcodes', right_on='rna_barcodes', how='left')
# %%

## filter out cells without RNA-cell type annotation
sub_bc_peak_count_adata = bc_peak_count_adata[~bc_peak_count_adata.obs['cell_type'].isna()].copy()
sub_entropy_bc_peak_count_adata = entropy_bc_peak_count_adata[~entropy_bc_peak_count_adata.obs['cell_type'].isna()].copy()
# %%

## On sub_bc_peak_count_adata, train PEAKVI model
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)
PEAKVI_LATENT_KEY = "X_peakvi"
PEAKVI_CLUSTERS_KEY = "clusters_peakvi"

#%%
scvi.model.PEAKVI.setup_anndata(sub_bc_peak_count_adata)
peak_model = scvi.model.PEAKVI(sub_bc_peak_count_adata)
peak_model.train()
# %%

# save latent representation back to peak_count_adata
latent = peak_model.get_latent_representation()
sub_bc_peak_count_adata.obsm[PEAKVI_LATENT_KEY] = latent
latent.shape
# save clustering results back to peak_count_adata
sc.pp.neighbors(sub_bc_peak_count_adata, use_rep=PEAKVI_LATENT_KEY)
sc.tl.umap(sub_bc_peak_count_adata, min_dist=0.2)
sc.tl.leiden(sub_bc_peak_count_adata, key_added=PEAKVI_CLUSTERS_KEY, resolution=0.2)

#%%
sc.pl.umap(sub_bc_peak_count_adata, color=[PEAKVI_CLUSTERS_KEY, 'cell_type'], title='latent space on atac peaks',)
# %%

## on sub_entropy_bc_peak_count_adata, train PEAKVI model
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)
PEAKVI_LATENT_KEY = "X_peakvi"
PEAKVI_CLUSTERS_KEY = "clusters_peakvi"

#%%
scvi.model.PEAKVI.setup_anndata(sub_entropy_bc_peak_count_adata)
entropy_peak_model = scvi.model.PEAKVI(sub_entropy_bc_peak_count_adata)
entropy_peak_model.train()
# %%

# save latent representation back to peak_count_adata
latent = entropy_peak_model.get_latent_representation()
sub_entropy_bc_peak_count_adata.obsm[PEAKVI_LATENT_KEY] = latent
latent.shape
# save clustering results back to peak_count_adata
sc.pp.neighbors(sub_entropy_bc_peak_count_adata, use_rep=PEAKVI_LATENT_KEY)
sc.tl.umap(sub_entropy_bc_peak_count_adata, min_dist=0.2)
sc.tl.leiden(sub_entropy_bc_peak_count_adata, key_added=PEAKVI_CLUSTERS_KEY, resolution=0.2)

#%%
sc.pl.umap(sub_entropy_bc_peak_count_adata, color=[PEAKVI_CLUSTERS_KEY, 'cell_type'], title='latent space on atac entropy peaks')
# %%
