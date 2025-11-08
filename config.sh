#!/bin/bash

export PYTHON_ENV='/home/syyang/python_virtuenv/atac_tobias3.13/bin/activate'
export PYTHON_ENV_PEAKVI='/home/syyang/python_virtuenv/mrvi_3.11/bin/activate'
export PYTHON_INVOKE_IN_R='/home/syyang/python_virtuenv/mrvi_3.11' # this python has to be built with '--enable-shared' option

export Genrich='/home/syyang/GitRepo/Genrich/./Genrich'

export Renv_Conda='seuratv5'

# make sure below tools available
# samtools
# genrich 
# bedtools 
# ripgrep (multi-thread grep) if filtering bam files on Barcodes that passed entropy threshold
# annotatePeaks.pl (from HOMER)




