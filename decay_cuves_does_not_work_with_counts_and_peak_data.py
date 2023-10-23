import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import curie as ci

def plot_decay_chain(parent_isotope, peak_data, units = 'd'):
    dc = ci.DecayChain(parent_isotope, A0=3e3, units=units)
    dc.counts = peak_data

    # dc.fit_A0()
    dc.plot(N_plot=2)







path =  '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/'
file_Ti01 = 'combined_peak_data_Ti01_copy.csv'
# path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/CJ010317_Ti01_18cm_30MeV/'

# file_Ti01 = 'CJ010317_Ti01_18cm_30MeV_peak_data_copy.csv'

file = file_Ti01

df = pd.read_csv(path+file,
        header=0,
        usecols=['isotope', 'start', 'stop', 'counts', 'unc_counts'])

df_Ti01_combined = df[['isotope', 'start', 'stop', 'counts', 'unc_counts']]

# df = df[['mean', '0', '1', '2', '3']]


print(df_Ti01_combined)



plot_decay_chain('48V', df_Ti01_combined)



# Measured counts: [start_time (d), stop_time (d), decays, unc_decays]
# Times relative to t=0 i.e. EoB time
# counts = {'225AC':[[5.0, 5.1, 6E5, 2E4],
#                           [6.0, 6.1, 7E5, 3E4]],
#                 '221FR':[5.5, 5.6, 6E5, 2E4]}
# def test_decay_curve(parent_isotope, counts, units = 'd'):
#     dc = ci.DecayChain(parent_isotope, units=units)
#     dc.counts=counts
#     dc.plot()

# test_decay_curve('225RA', counts)

