import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def fit_prod_rate(isotope, foil, path, file):
    R_estimated = 5
    t_irr_h = 0.33
    dc = ci.DecayChain(isotope, units='h', R=[[R_estimated, t_irr_h]])
    dc.get_counts(foil, '02/13/2017 14:27:00', path+file)

    isotopes, r, var_r = dc.fit_R()
    R = float(r[0])    
    var_R = float(var_r[0][0])
    dc.plot()
    return R, var_R


fit_prod_rate('56CO', 'Ni03', './', 'peak_data_onlyNi03_56Co.csv')
