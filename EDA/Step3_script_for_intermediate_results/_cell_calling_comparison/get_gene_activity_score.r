

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)

res_dir='/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_cell_calling_comparison'
fragment.path = paste0(res_dir, '/union_bc_fragments.tsv.gz')
counts_pic_file=paste0(res_dir, '/Entropy_cells_peaks/DownstreamReanalysis/bc_peak_pic_count.h5ad')
metadata_file=paste0(res_dir,'/Union_cell_3set.tsv')

library("anndata")
pic_h5 = read_h5ad(counts_pic_file) # cell-by-peak
counts_t = t(pic_h5$X)
colnames(counts_t) = pic_h5$obs[, 1]
rownames(counts_t) = pic_h5$var[, 1]

gene_activity_file = paste0(res_dir, '/Entropy_cells_peaks/DownstreamReanalysis/gene_activity_score.csv')


metadata <- read.csv(
  file = metadata_file,
  header = TRUE,
  row.names = 1,
  sep='\t'
)
rownames(metadata) = paste0(rownames(metadata), '-1')


chrom_assay <- CreateChromatinAssay(
  counts = counts_t,
  sep = c(":", "-"),
  fragments = fragment.path,
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


granges(pbmc)

# For hg38
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(pbmc) <- annotations

gene.activities <- GeneActivity(pbmc)
dense_gene.activities = as.matrix(gene.activities)


write.csv(dense_gene.activities, gene_activity_file, row.names = TRUE)


