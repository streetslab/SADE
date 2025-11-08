#%%
#### GLOBAL VARIABLES >>
precision = 1e-10 # for plot
#### GLOBAL VARIABLES << 

#%%



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot entropy distribution for cell barcodes and filtered fragments.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument('--chromosome', type=str, default='chr1', help='Chromosome used to calculate entropy.')
    parser.add_argument("--crbarcode_file", type=str, required=True, help="CellRanger barcode file, i.e., barcodes labeled to contain cells.")

    
    args = parser.parse_args()
    output_dir = args.output_dir
    crbarcode_file = args.crbarcode_file
    chromosome = args.chromosome
    # %%

    import pandas as pd
    import os
    import numpy as np

    import matplotlib.pyplot as plt
    # %%

    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy_df.tsv')
    bc_entropy = pd.read_csv(entropy_file, sep='\t', header=0, index_col=0)
    
    entropycutoff_file = os.path.join(output_dir, 'entropy_cutoff.csv')
    with open(entropycutoff_file, 'r') as f:
        line = f.readline()
        entropythreshold = float(line.strip().split(',')[1])
        
    bc_entropy['log10_entropy'] = np.log10(bc_entropy['Entropy'] + precision)
    bc_entropy['pass'] = bc_entropy['Entropy'] > entropythreshold
    

    # %%

    # >>>  compare entropy distribution between CR-labeled cells vs empty barcodes <<<
    cr_bc = pd.read_csv(crbarcode_file, sep='\t', header=None)

    cr_bc_entropy = bc_entropy.merge(cr_bc, left_index=True, right_on=0, how='right')
    cr_bc_entropy = cr_bc_entropy.set_index(0)

    cr_empty_entropy = bc_entropy[bc_entropy.index.isin(cr_bc[0]) == False]

    print("CellRanger labeled good vs bad barcodes stats:")
    print(f"    Number of good barcodes: {len(cr_bc_entropy)}")
    print(f"    Number of bad barcodes: {len(cr_empty_entropy)}")
    print(f"    Number of barcodes : {len(bc_entropy)} \n")

    with open(os.path.join(output_dir, 'entropy_stats.txt'), 'w') as f:
        f.write("Entropy statistics for CR_labeled cells:\n")
        cr_bc_entropy['Entropy'].describe().to_string(f)
        f.write('\n')
        f.write('entropythreshold: ' + str(entropythreshold) + '\n')
        f.write('\n')
        f.write("Entropy statistics for CR_empty barcodes:\n")
        cr_empty_entropy['Entropy'].describe().to_string(f)
        
    print("Entropy filtering good vs bad barcodes stats: ")
    print(f"    Number of good barcodes : {bc_entropy['pass'].sum()}")
    print(f"    Number of bad barcodes : {len(bc_entropy) - bc_entropy['pass'].sum()}")
    print(f"    Total number of barcodes : {len(bc_entropy)} \n")


    # %%
    # remove NaN values for plot
    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].violinplot(cr_bc_entropy['Entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[1].violinplot(cr_empty_entropy['Entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[0].set_ylabel('Entropy')
    ax[0].set_xlabel('CR cell barcodes')
    ax[1].set_xlabel('CR empty barcodes')
    ax[0].axhline(y=entropythreshold, color='r', linestyle='-', label=f'{entropythreshold:.4f}')
    ax[1].axhline(y=entropythreshold, color='r', linestyle='-', label=f'{entropythreshold:.4f}')
    fig.legend()
    fig_file = os.path.join(output_dir, 'figures', 'entropy_violin_EmptyCRcell.png')
    fig.savefig(fig_file)
    
    
    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].violinplot(cr_bc_entropy['log10_entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[1].violinplot(cr_empty_entropy['log10_entropy'].dropna(),  showmeans=True, showmedians=True, bw_method=0.1)
    ax[0].set_ylabel('log10 Entropy')
    ax[0].set_xlabel('CR cell barcodes')
    ax[1].set_xlabel('CR empty barcodes')
    ax[0].axhline(y=np.log10(entropythreshold ), color='r', linestyle='-', label=f'log10({entropythreshold:.4f})')
    ax[1].axhline(y=np.log10(entropythreshold ), color='r', linestyle='-', label=f'log10({entropythreshold:.4f})')
    fig.legend()
    fig_file = os.path.join(output_dir, 'figures', 'log10_entropy_violin_EmptyvsCell.png')
    fig.savefig(fig_file)

    # <<< Compare <<<


