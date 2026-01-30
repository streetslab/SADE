#!/bin/bash

export PYTHON_ENV='/home/syyang/python_virtuenv/atac_tobias3.13/bin/activate'
export PYTHON_ENV_PEAKVI='/home/syyang/python_virtuenv/mrvi_3.11/bin/activate'
export PYTHON_INVOKE_IN_R='/home/syyang/python_virtuenv/mrvi_3.11' # this python has to be built with '--enable-shared' option

export Genrich='/home/syyang/GitRepo/Genrich/./Genrich'

export Renv_Conda='seuratv5'

# Peak caller selection: "genrich" (requires BAM + samtools) or "macs" (uses fragments directly, no BAM needed)
export PEAK_CALLER='genrich'

# MACS version: "macs2" or "macs3" (only used when PEAK_CALLER=macs)
# macs3 recommended: active maintenance, HMMRATAC for ATAC-seq, scATAC-seq barcode support
export MACS_VERSION='macs3'

# make sure below tools available
# samtools (only needed for PEAK_CALLER=genrich)
# genrich (only needed for PEAK_CALLER=genrich)
# macs2/macs3 (only needed for PEAK_CALLER=macs; pip install MACS2 or MACS3)
# bedtools
# ripgrep (multi-thread grep) if filtering bam files on Barcodes that passed entropy threshold
# annotatePeaks.pl (from HOMER)
