
1. run sade on sample 
    bash sade.sh sample/fragments.tsv.gz -o VIB_10xultiome_2_WS3000F
    ## Able to run 
    ## /home/syyang/GitRepo/SADE/EDA/F1_demo_v.ipynb  

2. cd 'local_git_repo-manuscript_branch/EDA/reproduce_VIB_10xmultiome_2_analysisresult'
    - bash cp_scripts_to_resfolder.sh
        # Now you should have all the scripts in local_SADE-manuscript_branch/EDA/Step3_script_for_intermediate_results/* within VIB_10xultiome_2_WS3000F
    - cd VIB_10xultiome_2_WS3000F/_ArchR_TSS
    - bash prepare_frag_for_archr.sh 
        # Now you should have ArchR-TSSe filtering results, its filtered barcodes in the VIB_10xmultiome_2_WS3000F/_ArchR_TSS/ subfolder
    
    - cd VIB_10xultiome_2_WS3000F/_CR_FRIP
    - bash get_CR_barcodes.sh

    - run  /home/syyang/GitRepo/SADE/EDA/F2_prepare_data.ipynb
    ## running this is necessary to have 'Union_cell_3set_all_info.tsv' (BC that's not in union-bc are flagged but still kept) in the 'VIB_10xmultiome_2_WS3000F/_cell_calling_comparison'

    - cd VIB_10xultiome_2_WS3000F/_cell_calling_comparison
    - bash get_peak_sets.sh
    ## running this to have the Raw peaks, SADE-peaks, CR-peaks, and TSSe-peaks, 
    ##         as well as the 'Union_cell_3set_all_info.tsv' that only has the union-bc 
    ## Make sure this runs all the way through 
---- !!!!!
    - bash filter_frag.sh
    ## run this to have the union_bc_fragments.bed that has subset fragments corresponding to the union bc retained by any of the 3 methods
    ## -- here not run yet --
    - bash overlap_frag_peaks.sh 
    ## run this to have subset fragments that overlapped with each peak set.

    ## Able to run 
    ##  /home/syyang/GitRepo/SADE/EDA/F2_union_set_comparison.ipynb
    