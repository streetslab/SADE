#%%

Nsim = 1000  # number of simulated doublets
precision = 1e-10  # precision to avoid log(0)
        
        
        
#%% 
if __name__ == "__main__":
    
    #%% 
    import os
    import numpy as np
    import pandas as pd
    import argparse
    import pickle 
    #%%
    parser = argparse.ArgumentParser(description="Simulate doublets from single cell ATAC-seq data.")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to the output directory.")
    parser.add_argument("--chromosome", type=str, default="chr1", help="Chromosome to analyze.")


    output_dir = parser.parse_args().output_dir
    chromosome = parser.parse_args().chromosome
    
    #%%
    cluster_barcode_subdir = os.path.join(output_dir, 'cluster_barcodes')
    #  load insert frequency per barcode toå prep for doublet simulation
    insert_frequency_file = os.path.join(output_dir, f'{chromosome}_insert_frequency.pickle')
    with open(insert_frequency_file, 'rb') as file:
        insert_record = pickle.load(file)

    #%% 
    # Load barcodes for >= 2 populations
    cluster_barcode_file_list = os.listdir(cluster_barcode_subdir)
    
    if len(cluster_barcode_file_list) >= 2:
        ## provide barcode list for population 1    
        cluster1_barcode_file = os.path.join(cluster_barcode_subdir, cluster_barcode_file_list[0])
        cluster1_barcode_list = pd.read_csv(cluster1_barcode_file, header=None)[0].tolist()    
        ## provide barcode list for population 2
        cluster2_barcode_file = os.path.join(cluster_barcode_subdir, cluster_barcode_file_list[3])
        cluster2_barcode_list = pd.read_csv(cluster2_barcode_file, header=None)[0].tolist()
    else:
        raise ValueError("Need at least 2 populations to simulate doublets. \
                         TODO: implement random sampling on the same population doublets as negative control.")

    #%%
    ## Simulate doublets and calculate entropy 
    entropies_sim_doublets = []
    for i in range(Nsim):
        # randomly select one barcode from each population
        bc1 = np.random.choice(cluster1_barcode_list)
        bc2 = np.random.choice(cluster2_barcode_list)
        # simulate doublet insertions -- add insertion frequencies
        sim_insertions = insert_record[bc1].toarray() + insert_record[bc2].toarray()
        # calculate entropy 
        occur, freq = np.unique(sim_insertions, return_counts=True)
        freq = freq / np.sum(freq)
        entropy = -np.sum(freq * np.log2(freq + precision))
        #print(f'entropy is {entropy}')
        entropies_sim_doublets.append(entropy)
        
        
    #%%
    ## Simulate doublets from the same population as negative control
    entropies_sim_doublets_samepop = []
    for i in range(Nsim):
        # randomly select one barcode from each population
        bc1 = np.random.choice(cluster1_barcode_list)
        bc2 = np.random.choice(cluster1_barcode_list)
        # simulate doublet insertions -- add insertion frequencies
        sim_insertions = insert_record[bc1].toarray() + insert_record[bc2].toarray()
        # calculate entropy 
        occur, freq = np.unique(sim_insertions, return_counts=True)
        freq = freq / np.sum(freq)
        entropy = -np.sum(freq * np.log2(freq + precision))
        #print(f'entropy is {entropy}')
        entropies_sim_doublets_samepop.append(entropy)

    #%%
    # save simulated doublet entropies
    doublet_simulation_subdir = os.path.join(output_dir, 'doublet_simulation')
    if not os.path.exists(doublet_simulation_subdir):
        os.makedirs(doublet_simulation_subdir)
    sim_entropy_file = os.path.join(doublet_simulation_subdir, f'{chromosome}_simulated_doublet_entropies.csv')
    
    
    
    
    #%%
    import matplotlib.pyplot as plt
    figures_subdir = os.path.join(output_dir, 'figures')
    # Histogram of simulated doublet entropies
    fig, ax = plt.subplots(1, 1)
    ax.hist(np.log10(np.array(entropies_sim_doublets) + precision), bins=50, color='gray', alpha=0.7)
    ax.set_title(f'Simulated doublet entropies for {chromosome}')
    ax.set_xlabel('log10(Entropy)')
    ax.set_ylabel('Frequency')
    fig.savefig(os.path.join(figures_subdir, f'{chromosome}_simulated_doublet_entropies.png'))
    
    
    # Histogram of simulated doublet entropies from the same population
    fig, ax = plt.subplots(1, 1)
    ax.hist(np.log10(np.array(entropies_sim_doublets_samepop) + precision
                    ), bins=50, color='gray', alpha=0.7)
    ax.set_title(f'Simulated doublet entropies from the same population for {chromosome}')
    ax.set_xlabel('log10(Entropy)')
    ax.set_ylabel('Frequency')
    fig.savefig(os.path.join(figures_subdir, f'{chromosome}_simulated_doublet_entropies_samepop.png'))


    #%% 
    # Plot entropy distribution for all barcodes (cells + empty droplets) -- most of them are singlets
    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle')
    with open(entropy_file, 'rb') as f:
        entropy_dict = pickle.load(f)
        
    bc_entropies = list(entropy_dict.values())
    
    #%%
    fig, ax = plt.subplots(1, 1)
    ax.hist(np.log10(np.array(bc_entropies) + precision), bins=50, color='gray', alpha=0.7)
    ax.set_title(f'Barcode entropies for {chromosome}')
    ax.set_xlabel('log10(Entropy)')
    ax.set_ylabel('Frequency')


    ## plot them together
    fig, ax = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(8, 6))
    ax[0].hist(np.log10(np.array(bc_entropies) + precision), 
               bins=50, color='gray', alpha=0.7)
    ax[0].tick_params(labelbottom=True)
    ax[0].set_title(f'Barcode entropies for {chromosome}')
    ax[0].set_xlabel('log10(Entropy)')
    ax[0].set_ylabel('Frequency')
    
    ax[1].hist(np.log10(np.array(entropies_sim_doublets) + precision), bins=50, color='gray', alpha=0.7)
    ax[1].set_title(f'Simulated doublet entropies for {chromosome}')
    ax[1].set_ylabel('Frequency')
    ax[1].set_xlabel('log10(Entropy)')

    ax[2].hist(np.log10(np.array(entropies_sim_doublets_samepop) + precision), bins=50, color='gray', alpha=0.7)
    ax[2].set_title(f'Simulated doublet entropies from the same population for {chromosome}')
    ax[2].set_ylabel('Frequency')
    ax[2].set_xlabel('log10(Entropy)')
    
    fig.tight_layout()
    fig.savefig(os.path.join(figures_subdir, f'{chromosome}_simulated_doublet_entropies_comparison.png'))
    

    ## plot them together but freq is log scale for ref data 
    fig, ax = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(8, 6))
    ax[0].hist(np.log10(np.array(bc_entropies) + precision), 
               bins=50, color='gray', alpha=0.7, log=True)
    ax[0].tick_params(labelbottom=True)
    ax[0].set_title(f'Barcode entropies for {chromosome}')
    ax[0].set_xlabel('log10(Entropy)')
    ax[0].set_ylabel('Frequency')
    
    ax[1].hist(np.log10(np.array(entropies_sim_doublets) + precision), bins=50, color='gray', alpha=0.7)
    ax[1].set_title(f'Simulated doublet entropies for {chromosome}')
    ax[1].set_ylabel('Frequency')
    ax[1].set_xlabel('log10(Entropy)')

    ax[2].hist(np.log10(np.array(entropies_sim_doublets_samepop) + precision), bins=50, color='gray', alpha=0.7)
    ax[2].set_title(f'Simulated doublet entropies from the same population for {chromosome}')
    ax[2].set_ylabel('Frequency')
    ax[2].set_xlabel('log10(Entropy)')
    
    fig.tight_layout()
    fig.savefig(os.path.join(figures_subdir, f'{chromosome}_simulated_doublet_entropies_comparison2.png'))
# %%
