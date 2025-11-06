#%%
import numpy as np
import matplotlib.pyplot as plt



#%%
lam = 0.5

poi_sample = np.random.poisson(lam, size=1000)


#%%
val, freq = np.unique(poi_sample, return_counts=True)
freq_norm = freq / freq.sum()
entropy = -np.sum(freq_norm * np.log2(freq_norm))


fig, ax = plt.subplots()
ax.plot(val, freq_norm, marker='o')
ax.set_yscale('log')

# %%
