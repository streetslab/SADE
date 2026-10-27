# GRCh38-cCREs.bed is ENCODE SCREEN annotation regions (version: v4)
# https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed
wget "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed"
# Bed file description are: chr, start, end, accession_id1, accession_id2, cCRE_class
# accession_id1 and accession_id2 are documented in  https://github.com/weng-lab/ENCODE-cCREs/blob/master/Version-4/cCRE-Pipeline/Toolkit/make-region-accession.py
# these 2 columns have format  'EH38'+"mod_id"+"7 digits", where "mod_id" are one of the 3 letters ( 'D','E','F' )
# with modeDict = {"rDHS":"D","ccRE":"E", "cCRE":"E", "TF":"F"}
