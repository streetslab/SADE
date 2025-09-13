#### GLOBAL VARIABLES >>
EntropyThreshold = 0.000609
EntropyThreshold = 0.001
precision = 1e-10

#### GLOBAL VARIABLES << 



## see if barcodes filtering in CR would separate the entropy distributions within each sample




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot entropy distribution for cell barcodes and filtered fragments.")
    parser.add_argument("--res_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument("--entropy_file", type=str, required=True, help="pickle file of calculated entropies.")
    parser.add_argument("--crbarcode_file", type=str, required=True, help="CellRanger barcode file, i.e., barcodes labeled to contain cells.")
    parser.add_argument("--frag_file", type=str, required=True, help="Fragment file to filter on barcodes that passed entropy threshold.")
    parser.add_argument('--entropythreshold', type=float, default=EntropyThreshold, help='Entropy threshold for filtering barcodes.')
    
    
    args = parser.parse_args()
    res_dir = args.res_dir
    entropy_file = args.entropy_file
    crbarcode_file = args.crbarcode_file
    frag_file = args.frag_file
    entropythreshold = args.entropythreshold
    # %%

    import pandas as pd
    import pickle
    import os
    import numpy as np
    
    with open(entropy_file, 'rb') as f:
        entropy_dict = pickle.load(f)
    bc_entropy = pd.DataFrame.from_dict(entropy_dict, orient='index', columns=['entropy'])
    
    
    cr_bc = pd.read_csv(crbarcode_file, sep='\t', header=None)



    cr_bc_entropy = bc_entropy.merge(cr_bc, left_index=True, right_on=0, how='right')
    cr_bc_entropy = cr_bc_entropy.set_index(0)


    cr_empty_entropy = bc_entropy[bc_entropy.index.isin(cr_bc[0]) == False]
    
    print(f"Number of cell barcodes: {len(cr_bc_entropy)}")
    print(f"Number of empty barcodes: {len(cr_empty_entropy)}")
    print(f"Number of barcodes : {len(bc_entropy)}")
    
    with open(os.path.join(res_dir, 'entropy_stats.txt'), 'w') as f:
        f.write("Entropy statistics for CR_labeled cells:\n")
        cr_bc_entropy.describe().to_string(f)
        f.write('\n')
        f.write('entropythreshold: ' + str(entropythreshold) + '\n')
        f.write('\n')
        f.write("Entropy statistics for CR_empty barcodes:\n")
        cr_empty_entropy.describe().to_string(f)

    # %%
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].violinplot(cr_bc_entropy['entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
    ax[1].violinplot(cr_empty_entropy['entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
    ax[0].set_ylabel('Entropy')
    ax[0].set_xlabel('CR cell barcodes')
    ax[1].set_xlabel('CR empty barcodes')
    ax[0].axhline(y=entropythreshold, color='r', linestyle='-', label=f'{entropythreshold}')
    ax[1].axhline(y=entropythreshold, color='r', linestyle='-', label=f'{entropythreshold}')
    fig.legend()
    fig_file = os.path.join(res_dir, 'figures', 'entropy_violin_EmptyvsCell.png')
    fig.savefig(fig_file)
    
    
    cr_bc_entropy['log10_entropy'] = np.log10(cr_bc_entropy['entropy'] + precision)
    cr_empty_entropy['log10_entropy'] = np.log10(cr_empty_entropy['entropy'] + precision)
    fig, ax = plt.subplots(1, 2, sharey=True)
    ax[0].violinplot(cr_bc_entropy['log10_entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
    ax[1].violinplot(cr_empty_entropy['log10_entropy'],  showmeans=True, showmedians=True, bw_method=0.1)
    ax[0].set_ylabel('log10 Entropy')
    ax[0].set_xlabel('CR cell barcodes')
    ax[1].set_xlabel('CR empty barcodes')
    ax[0].axhline(y=np.log10(entropythreshold + precision), color='r', linestyle='-', label=f'log10({entropythreshold})')
    ax[1].axhline(y=np.log10(entropythreshold + precision), color='r', linestyle='-', label=f'log10({entropythreshold})')
    fig.legend()
    fig_file = os.path.join(res_dir, 'figures', 'log10_entropy_violin_EmptyvsCell.png')
    fig.savefig(fig_file)

    #%% 
    # now provide only BC that pass entropy threshold
    bc_entropy['pass'] = bc_entropy['entropy'] > entropythreshold
    
    print(f"Entropy threshold: {entropythreshold}")
    print(f"Number of barcodes that pass entropy threshold: {bc_entropy['pass'].sum()}")
    print(f"Number of barcodes that do not pass entropy threshold: {(bc_entropy['pass']==False).sum()}")
    
    
    with open(os.path.join(res_dir, 'stats_entropy_threshold.txt'), 'w') as f:
        f.write(f"Entropy threshold: {entropythreshold}\n")
        f.write(f"Number of barcodes that pass entropy threshold: {bc_entropy['pass'].sum()}\n")
        f.write(f"Number of barcodes that do not pass entropy threshold: {(bc_entropy['pass']==False).sum()}\n")
        f.write("\n")

    # save bc_entropy with pass column
    bc_pass_entropy_file = os.path.join(res_dir, 'bc_pass_entropy.tsv')
    bc_pass_entropy = bc_entropy.loc[bc_entropy['pass'] == True]
    bc_pass_entropy.to_csv(bc_pass_entropy_file, sep='\t')

    # %%


    import subprocess
    # filter fragment file based on cell barcodes that pass entropy threshold
    temp_bc_file = os.path.join(res_dir, "temp_bc.txt")
    filtered_frag_file = os.path.join(res_dir, "filtered_fragments.tsv") 
    subprocess.run(f"awk -v OFS=''  '{{print $1}}' {bc_pass_entropy_file} > {temp_bc_file}", shell=True)
    subprocess.run(f"grep -f {temp_bc_file} {frag_file} > {filtered_frag_file}", shell=True)
    subprocess.run(f"awk -v OFS='' -v prefix='CB:Z:' '{{print prefix, $1}}' {bc_pass_entropy_file} > {res_dir}/temp_bc_CBZ.txt", shell=True)


    frag_df = pd.read_csv(frag_file, sep='\t', header=None, comment='#')


    frag_df.columns = ['chrom', 'left_insert', 'right_insert', 'cb', 'support']
    frag_df['length'] = frag_df['right_insert'] - frag_df['left_insert']
    filter_frag_df = frag_df[frag_df['cb'].isin(bc_pass_entropy.index)]
    garbage_frag_df = frag_df[~ frag_df['cb'].isin(filter_frag_df['cb'])]
    # filter_frag_df = pd.read_csv(filtered_frag_file, sep='\t', header=None, comment='#')
    # filter_frag_df.columns = ['chrom', 'left_insert', 'right_insert', 'cb', 'support']
    # filter_frag_df['length'] = filter_frag_df['right_insert'] - filter_frag_df['left_insert']

    
    #%%

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 6), sharey=True)
    ax[0].hist(filter_frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[0].set_xlabel('Fragment length')
    ax[0].set_ylabel('Frequency')
    ax[0].set_title('Fragment length distribution for filtered fragments')
    ax[1].hist(frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[1].set_xlabel('Fragment length')
    ax[1].set_ylabel('Frequency')
    ax[1].set_title('Fragment length distribution for all fragments')
    ax[2].hist(garbage_frag_df['length'], bins=100, alpha=0.5, label='fl', color='green', range=(0,500))
    ax[2].set_xlabel('Fragment length')
    ax[2].set_ylabel('Frequency')
    ax[2].set_title('thrown-away fragments')
    fig_file = os.path.join(res_dir, 'figures', 'frag_df_length.png')
    fig.tight_layout()
    fig.savefig(fig_file)


