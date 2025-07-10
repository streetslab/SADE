#%% 
import anndata 
import pandas as pd
import numpy as np
import os 


import matplotlib.pyplot as plt
# %%
## load rna and atac barcode map 
multiome_dir = '/mnt/files/syy/adipose/multiom'

atac_rna_barcode_file = os.path.join(multiome_dir, 'atac_rna_barcodes_map.tsv')
atac_rna_barcode_map = pd.read_csv(atac_rna_barcode_file, sep='\t')

#%%
## load cell type annotation analyzed using RNA only 
rna_celltype_annotation_file = '/home/syyang/GitRepo/atac/Analysis_adipose/Henrique_Sean_RNAonly_celltype_multiome_samples.csv'

rna_celltype_annotation = pd.read_csv(rna_celltype_annotation_file, sep=',')

rna_celltype_annotation['rna_barcodes'] = rna_celltype_annotation['barcode'].apply(lambda x: x.split('-')[0])
#%%


# %%

Sample='MAPLE013_SC' 



sample_rna_celltype = rna_celltype_annotation[rna_celltype_annotation['sample'] == Sample].copy()
# merge by barcodes correspondance
sample_rna_merge = sample_rna_celltype.merge(atac_rna_barcode_map, left_on='rna_barcodes', right_on='rna_barcodes', how='left')


# %%
