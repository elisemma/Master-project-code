import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import curie as ci

def plot_decay_chain(parent_isotope, peak_data, units = 'd'):
    dc = ci.DecayChain(parent_isotope, A0=3E4, units=units)
    dc.counts = peak_data
    dc.fit_A0()
    print('A0: ',dc.A0)
    dc.plot()


file = 'peak_data.csv' 
# file = 'peak_data_start_time_as_mm_dd_osv.csv'
# file = 'peak_data_decays.csv'

df = pd.read_csv(file,
        header=0,
        usecols=['isotope', 'start', 'stop', 'counts', 'unc_counts'])
        # usecols=['isotope', 'start', 'stop', 'decays', 'unc_decays'])


print(df)

plot_decay_chain('48V', df)
