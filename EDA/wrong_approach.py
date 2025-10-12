## This was wrong -- but the plot looks great -- so might be useful to  see what's the reason
## Different distribution assumption??

# calcualte entropy without 0-insertion windows [Wrong approach -- made a mistake on calculating P0_open_region ]
for bc in tn5_insert_freq.keys():
    if bc not in entropy_df.index:
        continue
    else:
        cur, freq = np.unique(tn5_insert_freq[bc].toarray(), return_counts=True)
        if cur[0] == 0:
            freq_none = freq[1:]
            cur_none = cur[1:]
        else:
            freq_none = freq
            cur_none = cur


        total_windows = freq.sum() 
        

        non0_p = freq_none / np.sum(freq_none)
        non0_mean = np.sum(cur_none * non0_p)
        non0_2ndorder_moment = np.sum( cur_none**2 * non0_p)


        M1 = freq_none.sum()

        a = (non0_2ndorder_moment - non0_mean**2) / (non0_mean)   # made a mistake here, should swap non0_mean**2 and non0_mean


        if a  < precision:
            Nwindows_0inser_open_region = None
            entropy_final = None
        else:
            Nwindows_0inser_open_region = np.max((np.floor(M1 / a  - M1 ), 0))
            Nwindows_0inser_open_region = np.min((Nwindows_0inser_open_region, total_windows - M1))
            Nwindows_closed_region = int(Mtotal - M1 - Nwindows_0inser_open_region)

            total_ocur = np.concatenate((cur_none, [Nwindows_0inser_open_region, Nwindows_closed_region]))

            # calculate entropy
            posi_freq = total_ocur[total_ocur > 0]
            p_final = posi_freq / posi_freq.sum()

            entropy_final = -np.sum(p_final * np.log2(p_final))

        entropy_df.loc[bc, 'entropy_mixdist_mistake'] = entropy_final