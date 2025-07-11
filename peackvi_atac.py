#%%
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
import torch


def train_atac(adata_file, des_h5ad_file):
    
    scvi.settings.seed = 0
    print("Last run with scvi-tools version:", scvi.__version__)

    atac_adata = sc.read_h5ad(adata_file)
    # %%
    PEAKVI_LATENT_KEY = "X_peakvi"
    PEAKVI_CLUSTERS_KEY = "clusters_peakvi"

    #%%
    # filter regions with less than 3% of the cells
    # print("# regions before filtering:", bc_peak_count_adata.shape[-1])

    # sc.pp.filter_genes(bc_peak_count_adata, min_cells=10)
    # print("# regions after filtering:", bc_peak_count_adata.shape[-1])

    # train on atac_adata
    scvi.model.PEAKVI.setup_anndata(atac_adata)
    peak_model = scvi.model.PEAKVI(atac_adata)
    peak_model.train()
    # model_dir = os.path.join(save_dir.name, "peakvi_pbmc")
    # model.save(model_dir, overwrite=True)

    # save latent representation back adata
    latent = peak_model.get_latent_representation()
    atac_adata.obsm[PEAKVI_LATENT_KEY] = latent
    latent.shape
    # save clustering results back to atac_adata
    sc.pp.neighbors(atac_adata, use_rep=PEAKVI_LATENT_KEY)
    sc.tl.umap(atac_adata, min_dist=0.2)
    sc.tl.leiden(atac_adata, key_added=PEAKVI_CLUSTERS_KEY, resolution=0.2)

    # %%
    # save results
    atac_adata.write_h5ad(des_h5ad_file)





if __name__ == "__main__":
    import argparse
    import os
    arg = argparse.ArgumentParser(description="Run PEAKVI on ATAC data")
    arg.add_argument("--res_dir", type=str, required=True, help="Output directory for the results")
    args = arg.parse_args()

    output_dir = args.res_dir
    downstream_reanalysis_dir = os.path.join(output_dir, "DownstreamReanalysis")
    
    bc_peak_count_h5 = os.path.join(downstream_reanalysis_dir, "bc_peak_count.h5ad")
    bc_peak_peakVI_h5 = os.path.join(downstream_reanalysis_dir, "peakvi.h5ad")
    
    entropy_bc_peak_count_h5 = os.path.join(downstream_reanalysis_dir, 'entropy_bc_peak_count.h5ad')
    entropy_bc_peak_peakVI_h5 = os.path.join(downstream_reanalysis_dir, 'entropy_peakvi.h5ad')
    
    # train peakvi on bc_peak_count_h5
    train_atac(bc_peak_count_h5, bc_peak_peakVI_h5)
    
    # train peakvi on entropy_bc_peak_count_h5
    train_atac(entropy_bc_peak_count_h5, entropy_bc_peak_peakVI_h5)
    

    


