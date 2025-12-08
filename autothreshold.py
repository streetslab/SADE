import matplotlib.pyplot as plt
import numpy as np
import pickle


from scipy.interpolate import BSpline, make_splrep


def fit_spline_and_find_cutoff(x: np.ndarray, y: np.ndarray, k: int=2, spline_s: int=67, spline_k: int=3, limit: int=20000) -> tuple[float, int, BSpline]:
    """
    Takes the ranks (x) and the corresponding log 10 entropy values (y) and returns estimated entropy cutoff, the rank cutoff, and the scipy BSpline object used to fit the curve.
    The function fits a smoothing spline to the curve, then computes the first and second derivatives on the spline. then iterating on increasing x values, finds regions where the left bound of the region is where second derivative transitions from greater than cutoff to below cutoff, and then the right bound is where second derivative goes from below cutoff to above. then takes the kth such region by ranks, and computes the greatest dropoff in the region (min second derivative). returns this point as cutoff 

    x: np.ndarray, the ranks (1 is highest)
    y: np.ndarray, the log_10 entropy values
    k: int, kth region where second derivative transitions from > cutoff to less than cutoff, then from less cutoff to greater than
    cutoff: float, close to zero, look at the derivative curves because if at 0 then too sensitive to tiny changes
    spline_s: int, spline smoothing conditions, see scipy.make_splrep
    spline_k: int, spline degree, see scipy.make_splrep

    returns: tuple[float, int, Bspline]
    first value is entropy cutoff, second value is rank cutoff, third is the spline fit to the curve and used to find cutoffs
    """
    # only consider up to limit barcodes for spline fitting (as even the highest throughput datasets have < 20k cells)
    x = x[:limit]
    y = y[:limit]
    cs = make_splrep(x, y, k=spline_k, s=spline_s)
    
    first_deriv = cs.derivative(1)(x)
    second_deriv = cs.derivative(2)(x)
    second_deriv = second_deriv / np.abs(second_deriv).max()  # normalize second derivative to better find inflection regions
    
    left = 0
    right = 0
    flag = False
    count = 0 

    inflection_regions = {}
    
    for i in x: 
        
        d_2 = second_deriv[i]
        if (not flag) and d_2 < 0:
            flag = True
            left = i
        elif flag and d_2 > 0:
            flag = False
            right = i
            
            inflection_regions[count] = (left, right) # store inflection region boundaries
            count += 1
        
        if count == k: 
            break

    # find the global minimum first derivative among the k inflection regions found
    regions_rank = []
    for region_index, (region) in inflection_regions.items():
        left, right = region
        if left > right:
            right = x[-1]
        rank = np.argmin(first_deriv[left:right+1]) + left
        regions_rank.append(rank)

    rank = regions_rank[np.argmin(first_deriv[regions_rank])]
    
    return rank, y[rank], cs

def get_x_y_from_pickle_helper(pickle_path: str):
    with open(pickle_path, 'rb') as f:
        entropies = pickle.load(f)

    sorted_entropies = dict(sorted(entropies.items(), key=lambda item: item[1], reverse=True))
    barcodes = sorted_entropies.keys()
    values = sorted_entropies.values()
    rank = range(len(barcodes))
    return np.array(list(rank)), np.array(np.log10(list(values)))

def make_plots(pickle_path, save_path):
    x, y = get_x_y_from_pickle_helper(pickle_path)
    rank_cutoff, entr_cutoff, cs = fit_spline_and_find_cutoff(x, y)

    spl_y = cs(x)
    deriv = {}
    for i in range(1, 3):
        deriv[i] = cs.derivative(i)(x)

    deriv[2] = deriv[2] / np.abs(deriv[2]).max()  # normalize second derivative for better visualization    
    
    fig, ax = plt.subplots(3, 1, figsize=(10,10), dpi=300)
    ax[0].plot(x, y, 'o', label='data')
    ax[0].plot(x, spl_y, label='spline')
    ax[1].plot(x, deriv[1], label='first derivative')
    ax[2].plot(x, deriv[2], label='second derivative')

    d1_limit = np.abs(deriv[1]).max()
    ax[1].set_ylim(-d1_limit*1.1, d1_limit*1.1)

    for i in range(3):
        ax[i].axvline(rank_cutoff, color='pink', linestyle='--', label='rank_cutoff')
        if i == 0:
            ax[i].axhline(entr_cutoff, color='pink', linestyle='--', label='entr_cutoff')

    ax[0].set_xlabel('rank')
    ax[0].set_ylabel('log 10 entropy')
    ax[0].set_title('fit spline')
    ax[1].set_title('first derivative')
    ax[2].set_title('second derivative')

    plt.tight_layout()
    fig.savefig(save_path)
