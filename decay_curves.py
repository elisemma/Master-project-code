import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import curie as ci

def plot_dacay_chain(parent_isotope, spectrum_list, peak_data, units = 'd'):

    
    dc = ci.DecayChain(parent_isotope, A0 = 3e3, units=units)
    # dc.plot()
    # dc.get_counts(spectrum_list, peak_data=peak_data, EoB='01/01/2016 08:39:08')

    dc.get_counts(spectrum_list, peak_data=peak_data, EoB='02/13/2017 08:39:08') #Not do this!!! and the time is wrong

    dc.fit_A0()
    dc.plot(N_plot=2)




path_spectrum = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/30MeV/'
path_peak_data_CJ = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/CJ010317_Ti01_18cm_30MeV/'
path_peak_data_CS = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/CS060317_Ti01_18cm_30MeV/'
spectrum_list = [path_spectrum+'CJ010317_Ti01_18cm_30MeV.Spe', path_spectrum+'CS060317_Ti01_18cm_30MeV.Spe']
# spectrum_list = path_spectrum+'CJ010317_Ti01_18cm_30MeV.Spe'
# peak_data = path_peak_data_CJ+'CJ010317_Ti01_18cm_30MeV_peak_data.csv'
peak_data = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/combined_peak_data_Ti01.csv'

plot_dacay_chain('48V', spectrum_list, peak_data)

















# path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/CJ010317_Ti01_18cm_30MeV/'
# file = 'CJ010317_Ti01_18cm_30MeV_peak_data.csv'

# df = pd.read_csv(path+file,
#         header=0,
#         usecols=['decays', 'unc_decays', 'start_time', 'live_time'])
# # printing dataframe
# print(df.loc('58V'))

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

