#%%
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
import torch


def train_atac(adata_file, des_h5ad_file, des_model_dir):
    
    scvi.settings.seed = 0
    print("Last run with scvi-tools version:", scvi.__version__)

    atac_adata = sc.read_h5ad(adata_file)
    # %%
    PEAKVI_LATENT_KEY = "X_peakvi"
    PEAKVI_CLUSTERS_KEY = "clusters_peakvi"

    #%%
    # filter regions with less than 3% of the cells
    # print("# regions before filtering:", peak_count_adata.shape[-1])

    # sc.pp.filter_genes(peak_count_adata, min_cells=10)
    # print("# regions after filtering:", peak_count_adata.shape[-1])

    # Let's Use All the identified peaks for comarison (since entropy approach already has more peaks)
    # train on atac_adata
    scvi.model.PEAKVI.setup_anndata(atac_adata)
    peak_model = scvi.model.PEAKVI(atac_adata)
    peak_model.train()
    # save model
    peak_model.save(des_model_dir, overwrite=True)

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
    standard_peak_peakVI_h5 = os.path.join(downstream_reanalysis_dir, "standard_peakvi.h5ad")
    standard_peakvi_model_dir = os.path.join(downstream_reanalysis_dir, "standard_peak_peakvi_model")
    
    entropy_bc_peak_count_h5 = os.path.join(downstream_reanalysis_dir, 'entropy_bc_peak_count.h5ad')
    entropy_peak_peakVI_h5 = os.path.join(downstream_reanalysis_dir, 'entropy_peakvi.h5ad')
    entropy_peakvi_model_dir = os.path.join(downstream_reanalysis_dir, "entropy_peak_peakvi_model")

    # train peakvi on standard_peak_count_h5
    train_atac(bc_peak_count_h5, standard_peak_peakVI_h5, standard_peakvi_model_dir)
    
    # train peakvi on entropy_peak_count_h5
    train_atac(entropy_bc_peak_count_h5, entropy_peak_peakVI_h5, entropy_peakvi_model_dir)
    


