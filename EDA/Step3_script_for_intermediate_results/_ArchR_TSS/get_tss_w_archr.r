
library(ArchR)
packageVersion("ArchR") #[1] ‘1.0.3’
library(parallel)
set.seed(1)

inputFiles=c('/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_ArchR_TSS/fragments.tsv.gz') # copied to the current working directory from the original fragments.tsv.gz file
names(inputFiles)='VIB_10xmultiome_2'

addArchRGenome("hg38")
addArchRThreads(threads = 16) 

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


## Save ArchR results to dataframe
metadata=readRDS("/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_ArchR_TSS/QualityControl/VIB_10xmultiome_2/VIB_10xmultiome_2-Pre-Filter-Metadata.rds")
write.table(metadata, sep='\t', '/mnt/hdd_bob/syy/adipose/atac/res/VIB_10xmultiome_2_WS3000F/_ArchR_TSS/QualityControl/VIB_10xmultiome_2_Metadata.tsv', row.names=FALSE)
