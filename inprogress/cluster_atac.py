#%%
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
import torch
# %%

bc_peak_count_h5 = '/mnt/hdd_bob/syy/adipose/atac/res/NK_ATAC_MAPLE013_SC/DownstreamReanalysis/bc_peak_count.h5ad'
bc_peak_peakVI_h5 = '/mnt/hdd_bob/syy/adipose/atac/res/NK_ATAC_MAPLE013_SC/DownstreamReanalysis/peakvi.h5ad'

# %%
bc_peak_count_adata = sc.read_h5ad(bc_peak_count_h5)

# %%

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)


# %%
PEAKVI_LATENT_KEY = "X_peakvi"
PEAKVI_CLUSTERS_KEY = "clusters_peakvi"

#%%
# filter regions with less than 3% of the cells
print("# regions before filtering:", bc_peak_count_adata.shape[-1])

sc.pp.filter_genes(bc_peak_count_adata, min_cells=10)
print("# regions after filtering:", bc_peak_count_adata.shape[-1])

# train PEAKVI model on bc_peak_count_adata
scvi.model.PEAKVI.setup_anndata(bc_peak_count_adata)
peak_model = scvi.model.PEAKVI(bc_peak_count_adata)
peak_model.train()

# model_dir = os.path.join(save_dir.name, "peakvi_pbmc")
# model.save(model_dir, overwrite=True)

# %%

# save latent representation back to peak_count_adata
latent = peak_model.get_latent_representation()
bc_peak_count_adata.obsm[PEAKVI_LATENT_KEY] = latent
latent.shape
# save clustering results back to peak_count_adata
sc.pp.neighbors(bc_peak_count_adata, use_rep=PEAKVI_LATENT_KEY)
sc.tl.umap(bc_peak_count_adata, min_dist=0.2)
sc.tl.leiden(bc_peak_count_adata, key_added=PEAKVI_CLUSTERS_KEY, resolution=0.2)

#%%
sc.pl.umap(bc_peak_count_adata, color=PEAKVI_CLUSTERS_KEY)
# %%
# save results
bc_peak_count_adata.write_h5ad(bc_peak_peakVI_h5)

# %%
