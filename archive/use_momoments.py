##


def calculate_moment(tn5_inser_array:np.ndarray, return_freq:bool=False, print_debug:bool=False):
    cur, freq = np.unique(tn5_inser_array, return_counts=True)
    if print_debug:
        print(cur, freq)
    if cur[0] < precision:
        freq_none = freq[1:]
        cur_none = cur[1:]
    else:
        freq_none = freq
        cur_none = cur

    Mtotal = freq.sum() 
    M1 = freq_none.sum()
    
    non0_p = freq_none / np.sum(freq_none)
    non0_mean = np.sum(cur_none * non0_p)
    non0_2ndorder_moment = np.sum( cur_none**2 * non0_p)
    
    P0 = 1 - (non0_2ndorder_moment - non0_mean) / non0_mean**2

    if return_freq:
        return Mtotal, M1, non0_mean, non0_2ndorder_moment, P0, cur_none, freq_none
    
    return Mtotal, M1, non0_mean, non0_2ndorder_moment, P0


def calculate_moment_fix(tn5_inser_array:np.ndarray, return_freq:bool=False, print_debug:bool=False):
    cur, freq = np.unique(tn5_inser_array, return_counts=True)
    if print_debug:
        print(cur, freq)
    if cur[0] < precision:
        freq_none = freq[1:]
        cur_none = cur[1:]
    else:
        freq_none = freq
        cur_none = cur

    cur_none = (cur_none + 1) // 2 # adjust for paired-end sequencing
    
    Mtotal = freq.sum() 
    M1 = freq_none.sum()
    
    non0_p = freq_none / np.sum(freq_none)
    non0_mean = np.sum(cur_none * non0_p)
    non0_2ndorder_moment = np.sum( cur_none**2 * non0_p)
    
    P0 = 1 - (non0_2ndorder_moment - non0_mean) / non0_mean**2

    if return_freq:
        return Mtotal, M1, non0_mean, non0_2ndorder_moment, P0, cur_none, freq_none

    return Mtotal, M1, non0_mean, non0_2ndorder_moment, P0



# calcualte entropy without 0-insertion windows 
def calculate_entropy(tn5_insertion_array:np.ndarray):
    # tn5_insert_freq[bc].toarray()
    Mtotal, M1, non0_mean, non0_2ndorder_moment, p0_open_region, cur_none, freq_none = calculate_moment(tn5_insertion_array, return_freq=True)

    if p0_open_region < 0 or p0_open_region >= 1:
        Nwindows_0inser_open_region = None
    else:
        Nwindows_0inser_open_region = np.floor(p0_open_region / (1 - p0_open_region) * M1)

        Nwindows_0inser_open_region = np.max((Nwindows_0inser_open_region, Mtotal - M1))

    if Nwindows_0inser_open_region is None:
        Entropy_final = None
    else:
        Nwindows_closed_region = int(Mtotal - M1 - Nwindows_0inser_open_region)
        total_freq = np.concatenate(([Nwindows_closed_region, Nwindows_0inser_open_region], freq_none))

        # calculate entropy
        posi_freq = total_freq[total_freq > 0]
        p_final = posi_freq / posi_freq.sum()
        Entropy_final = -np.sum(p_final * np.log2(p_final))
    return Entropy_final


## Calcualte multiple entropies: entropy of mixture distribution, entropy of open region, entropy of final distribution
def calcualte_multiple_entropies(tn5_insert_array:np.ndarray, return_all_paras:bool=False):
    Mtotal, M1, non0_mean, non0_2ndorder_moment, p0_open_region, cur_none, freq_none = calculate_moment_fix(tn5_insert_array, return_freq=True)

    if p0_open_region <= 0 or p0_open_region >= 1:
        Nwindows_0inser_open_region = None
        Entropy_mixturedist = None
        Entropy_open_region = None
        Entropy_final = None
    else:
        Nwindows_0inser_open_region = np.floor(M1 / ( 1- p0_open_region) * p0_open_region )
        Nwindows_0inser_open_region = np.min((Nwindows_0inser_open_region, Mtotal - M1))
        Nwindows_closed_region = int(Mtotal - M1 - Nwindows_0inser_open_region)

        total_freq = np.concatenate(([Nwindows_closed_region, Nwindows_0inser_open_region], freq_none))

        # calculate entropy of mix states
        p_stats = np.array([Nwindows_closed_region, Mtotal - Nwindows_closed_region]) / Mtotal
        Entropy_state = -np.sum(p_stats * np.log2(p_stats))

        open_region_freq = total_freq[1:]
        open_region_freq_non0 = open_region_freq[open_region_freq > 0]
        open_region_p = open_region_freq_non0 / open_region_freq_non0.sum()

        # calculate entropy as a mixture distribution 
        Entropy_open_region = -np.sum(open_region_p * np.log2(open_region_p))
        Entropy_mixturedist = Entropy_state + (1 - Nwindows_closed_region / Mtotal) * Entropy_open_region

        # calculate entropy of final distribution -- Simply breaking down 0-freq to from open and closed regions
        posi_freq = total_freq[total_freq > 0]
        p_final = posi_freq / posi_freq.sum()
        Entropy_final = -np.sum(p_final * np.log2(p_final))

    if not return_all_paras:
        return Entropy_mixturedist, Entropy_open_region, Entropy_final
    
    other_paras = (non0_mean, non0_2ndorder_moment, p0_open_region)
    return Entropy_mixturedist, Entropy_open_region, Entropy_final, other_paras


def break_down_freq(tn5_insert_array:np.ndarray):
    cur, freq = np.unique(tn5_insert_array, return_counts=True)
    if cur[0] == 0:
        freq_none = freq[1:]
        cur_none = cur[1:]
    else:
        freq_none = freq
        cur_none = cur

    Mtotal = freq.sum() 
    M1 = freq_none.sum()
    
    non0_p = freq_none / np.sum(freq_none)
    non0_mean = np.sum(cur_none * non0_p)
    non0_2ndorder_moment = np.sum( cur_none**2 * non0_p)


    p0_open_region = 1 - (non0_2ndorder_moment - non0_mean) / non0_mean**2

    if p0_open_region < 0 or p0_open_region >= 1:
        Nwindows_0inser_open_region = None
    else:
        Nwindows_0inser_open_region = np.floor(p0_open_region / (1 - p0_open_region) * M1)
        Nwindows_0inser_open_region = np.max((Nwindows_0inser_open_region, Mtotal - M1))

    if Nwindows_0inser_open_region is None:
        total_freq = np.concatenate(([None, None], freq_none))
    else:
        Nwindows_closed_region = int(Mtotal - M1 - Nwindows_0inser_open_region)
        total_freq = np.concatenate(([Nwindows_closed_region, Nwindows_0inser_open_region], freq_none))
    
    total_cur = np.concatenate(([-1, 0], cur_none))

    return total_cur, total_freq

################