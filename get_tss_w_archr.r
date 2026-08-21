
suppressMessages(library("optparse"))
library(ArchR)
packageVersion("ArchR") #[1] ‘1.0.3’
library(parallel)
set.seed(1)


option_list <- list( 
    make_option(c("-f", "--fragments_file_in_gz"), type="character", default=FALSE,
              help="Fragments file in gzip format."),
    make_option(c("-o", "--res_dir"), type="character", default=FALSE,
              help="Directory to results."), 
    make_option(c("-n", "--sample_name"), type="character", default=FALSE,
              help="Sample name"),
    make_option(c("-g", "--genome_name"), type="character", default=FALSE,
              help="Genome name, e.g., mm10 or hg38."
    ))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fragments_file = opt$fragments_file_in_gz
output_dir = opt$res_dir
sample_name = opt$sample_name
genome_name = opt$genome_name

print(paste0("Fragments file: ", fragments_file))
print(paste0("Output directory: ", output_dir))
print(paste0("Sample name: ", sample_name))
print(paste0("Genome name: ", genome_name))

inputFiles=c(fragments_file)
names(inputFiles)=sample_name

addArchRGenome(genome_name)
addArchRThreads(threads = 32) 

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, # default from ArchR is 4
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


## Save ArchR results to dataframe
# > ArchR default save its results in such sub-dir structure as .rds file. So need to load it and save as tsv
archr_res_file=file.path(output_dir, "_ArchR_TSS/QualityControl", sample_name, paste0(sample_name, "-Pre-Filter-Metadata.rds"))
archr_tsv_filename=file.path(output_dir, "_ArchR_TSS/QualityControl", paste0(sample_name, "_Metadata.tsv"))

metadata = readRDS(archr_res_file)
write.table(metadata, sep='\t', archr_tsv_filename,  row.names=FALSE)
