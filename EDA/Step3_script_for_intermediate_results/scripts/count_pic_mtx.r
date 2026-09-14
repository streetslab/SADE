
#%%
suppressMessages(library("GenomicRanges"))
suppressMessages(library("Matrix"))
suppressMessages(library("PICsnATAC"))
suppressMessages(library("optparse"))
suppressMessages(library("anndata"))
suppressMessages(library(reticulate))
use_virtualenv("/home/syyang/python_virtuenv/mrvi_3.11", required = TRUE) # change to your own python virtual environment path that has mrvi/scvi package installed

#%%
option_list <- list( 
    make_option(c("-b", "--barcode_file"), type="character", default=FALSE,
              help="Barcode file."),
    make_option(c("-o", "--res_dir"), type="character", default=FALSE,
              help="Directory to peak-set"),
    make_option(c("-f", "--frag_file"), type="character", default=FALSE,
              help="fragment file tsv.gz")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


barcode_file = opt$barcode_file
filtered_frag_gz_file = opt$frag_file
print(paste("frag file", filtered_frag_gz_file))

output_dir = opt$res_dir
peak_file = file.path(output_dir, "peaks.bed")

DownstreamReanalysis_dir = file.path(output_dir, "DownstreamReanalysis")
if (!file.exists(DownstreamReanalysis_dir)) {
	  dir.create(DownstreamReanalysis_dir)
}
bc_peak_count_h5 = file.path(DownstreamReanalysis_dir, "bc_peak_pic_count.h5ad")

#%%
# Load barcodes 
barcodes_df = read.csv(barcode_file, sep='\t', header=FALSE, col.names='barcodes')
print(paste("Number of barcodes:", nrow(barcodes_df), ' barcodes will be used for PIC counting.'))
barcodes = barcodes_df$barcodes


# Load peak sets
peak_df = read.csv(peak_file, header=FALSE, sep='\t')
if (ncol(peak_df) < 3) {
    stop("The peak file must have at least 3 columns.")
}
colnames(peak_df) <- c("seqname", "start", "end") # 'seqname' is chromosome name
peak_df = peak_df[, 1:3]
peak_sets <- GenomicRanges::makeGRangesFromDataFrame(peak_df)


#%%
# Get the PIC matrix [cells x peaks]
pic_matrix <- PIC_counting(cells=barcodes, 
    fragment_tsv_gz_file_location=filtered_frag_gz_file, 
    peak_sets=peak_sets)

pic_matrix_adata = AnnData(X=t(pic_matrix), 
    obs = data.frame(barcodes=colnames(pic_matrix)),
    var = data.frame(peak=rownames(pic_matrix)))

pic_matrix_adata$write_h5ad(bc_peak_count_h5)

