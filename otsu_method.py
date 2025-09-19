#%%
import numpy as np 
import pickle 
import os 
import pandas as pd 
precision = 1e-10  # to avoid log(entropy==0) issues

#%%



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Otsu's method to determine entropy threshold")
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the output threshold file')
    parser.add_argument('--chromosome', type=str, default='chr1', help='Chromosome to analyze')
    args = parser.parse_args()


    output_dir = args.output_dir
    chromosome = args.chromosome
    
    # Load entropy data
    entropy_file = os.path.join(output_dir, f'{chromosome}_barcode_entropy.pickle') 
    with open(entropy_file, 'rb') as f:
        entropy_dic = pickle.load(f)

    entropy_values = np.array(list(entropy_dic.values())) 
    log10entropy_values = np.log10(entropy_values + precision)  # add precision to avoid log(0)

    # Compute histogram
    hist, bin_edges = np.histogram(log10entropy_values, bins=265, range=(-2.5, np.max(log10entropy_values)))  

    print('range is ', (-2.5, np.max(log10entropy_values)))
    # use log-frequency to present the histogram 
    
    # Calculate probabilities
    hist = hist.astype(float)
    # prob = hist / np.sum(hist)
    log_hist = np.log10(hist + precision)
    prob = log_hist / np.sum(log_hist)

    # Cumulative sums
    cumulative_prob = np.cumsum(prob)
    cumulative_mean1 = np.cumsum(prob * bin_edges[:-1])

    # Global mean
    global_mean = cumulative_mean1[-1]

    # Between-class variance
    sigma_b_squared = cumulative_prob / ((1 - cumulative_prob) + precision ) * (global_mean  - cumulative_mean1) ** 2
    # mean2 = (global_mean - cumulative_mean1 * cumulative_prob) / ( 1- cumulative_prob + precision)
    # sigma_b_squared = cumulative_prob * (1 - cumulative_prob) * (mean2 - cumulative_mean1) ** 2

    # Find the threshold that maximizes the between-class variance
    optimal_threshold_index = np.argmax(sigma_b_squared)
    optimal_threshold = bin_edges[optimal_threshold_index]

    print(f"Optimal log10-Entropy Threshold: {optimal_threshold}")

    # Save the threshold to a file