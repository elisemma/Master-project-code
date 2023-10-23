import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import curie as ci

def plot_decay_chain(parent_isotope, peak_data, units = 'd'):
    dc = ci.DecayChain(parent_isotope, R=[[5,1], [8,2]], units=units)
    # dc.counts = peak_data
    dc.counts = {'48V':[[5.0, 5.1, 6E5, 2E4], [6.0, 6.1, 7E5, 3E4]]}

    # decay = dc.decays('48V', 5.0, 5.1, 'd')
    # print('_______________________')
    # print(decay)
    # print('_______________________')

    dc.fit_R()
    print('R: ',dc.R)
    dc.plot()







path =  '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/'
file_Ti01 = 'temp.csv'
# path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/CJ010317_Ti01_18cm_30MeV/'

# file_Ti01 = 'CJ010317_Ti01_18cm_30MeV_peak_data_copy.csv'

file = file_Ti01

df = pd.read_csv(path+file,
        header=0,
        # usecols=['decays', 'unc_decays', 'live_time', 'real_time'])
        # usecols=['decays', 'unc_decays', 'start_time', 'live_time'])
        # usecols=['start_time', 'live_time', 'decays', 'unc_decays'])
        usecols=['isotope', 'start', 'stop', 'counts', 'unc_counts'])

print(df)

# df_Ti01_combined = df[['isotope', 'start', 'stop', 'counts', 'unc_counts']]
# print(df_Ti01_combined)

# df = df[['mean', '0', '1', '2', '3']]




# plot_decay_chain('48V', df_Ti01_combined)
plot_decay_chain('48V', df)



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

