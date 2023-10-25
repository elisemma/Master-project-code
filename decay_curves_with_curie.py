import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def fit_prod_rate(isotope, foil, path, file):
    R_estimated = 5
    t_irr_h = 0.33
    dc = ci.DecayChain(isotope, units='h', R=[[R_estimated, t_irr_h]])
    dc.get_counts(foil, '02/13/2017 14:27:00', path+file)

    isotopes, R, var_R = dc.fit_R()
    dc.plot(titel = foil)




isotope = '48V'
foil = 'Ti01'
path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/'
file = 'combined_peak_data_Ti01.csv'

fit_prod_rate(isotope, foil, path, file)