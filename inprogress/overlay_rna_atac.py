#%%

import scanpy as sc
import numpy as np
import pandas as pd
import os 

#%%

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Overlay RNA cell type annotation on ATAC peak latent space")
    parser.add_argument('--res_dir', type=str, required=True, help='Directory containing the results')
    args = parser.parse_args()
    
    output_dir = args.res_dir


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
    rna_celltype_annotation_file = '/home/syyang/GitRepo/atac/Analysis_adipose/Henrique_Sean_RNAonly_celltype_multiome_samples.csv'
    rna_celltype_annotation = pd.read_csv(rna_celltype_annotation_file, sep=',', index_col=0)
    rna_celltype_annotation['Sample'] = rna_celltype_annotation['sample'].apply(lambda x: x.upper()) # make it sample name upper case
    rna_celltype_annotation['rna_barcodes'] = rna_celltype_annotation['barcode'].apply(lambda x: x.split('-')[0])
    # %%

    sample_rna_celltype = rna_celltype_annotation[rna_celltype_annotation['Sample'] == rna_sample_name].copy()

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
    peakvi_fig_file = os.path.join(figure_dir, 'umap_peakvi.png')
    sc.pl.umap(peakvi_adata, color=['cell_type', 'clusters_peakvi'], ncols=1,  title='RNA-cell-type on atac peaks latent space', \
            size=30, na_color='None', save='umap_peakvi.png')
    entropy_peakvi_fig_file = os.path.join(figure_dir, 'umap_entropy_peakvi.png')
    sc.pl.umap(entropy_peakvi_adata, color=['cell_type', 'clusters_peakvi'], ncols=1,  title='RNA-cell-type on atac peaks latent space',  \
            size=30, na_color='None', save='umap_entropy_peakvi.png')

# %%
