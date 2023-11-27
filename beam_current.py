from statistics import variance
import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 
import glob
from scipy.interpolate import CubicSpline
import seaborn as sns


def convert_production_rate_to_A0(decay_const, prod_rate, prod_rate_unc, t_irr, t_irr_unc): 
    A0 = decay_const*prod_rate*t_irr
    # A0_unc = 
    return A0#, A0_unc




def convert_A0_to_beam_current_w_unc(A0, areal_dens, Mm, decay_const, xs_E_mon, A0_unc, areal_dens_unc_percent, Mm_unc, xs_E_mon_unc):
    N_A = 6.0221408e+23
    t_irr = 39240 #[s]
    t_irr_unc = 5 #[s]
    N_T = float(areal_dens)*N_A/Mm*10 #[nuclei/m^2] når den tar inn areal_dens i mg/cm^2

    beam_current = A0/(N_T*xs_E_mon*(1-np.exp(-decay_const*t_irr))) #[d/s]
    beam_current_in_A = beam_current*1.60217634e-19 #[A]
    beam_current_in_nA = beam_current_in_A*1e9 #[nA]

    areal_dens_unc = areal_dens*10*areal_dens_unc_percent/100
    N_T_unc = N_T*np.sqrt((areal_dens_unc/areal_dens)**2 + (Mm_unc/Mm)**2) #[nuclei/cm^2]

    dfdx_list = [] #Jacobian
    unc_list = []

    dA0 = A0*1e-8
    dfdA0 = (A0+dA0/(N_T*xs_E_mon*(1-np.exp(-decay_const*t_irr))) - A0-dA0/(N_T*xs_E_mon*(1-np.exp(-decay_const*t_irr))))/dA0
    dfdx_list.append(dfdA0)
    unc_list.append(A0_unc)

    dN_T = N_T*1e-8
    dfdN_T = (A0/((N_T+dN_T)*xs_E_mon*(1-np.exp(-decay_const*t_irr))) - A0/((N_T-dN_T)*xs_E_mon*(1-np.exp(-decay_const*t_irr))))/dN_T
    dfdx_list.append(dfdN_T)
    unc_list.append(N_T_unc)

    dxs = xs_E_mon*1e-8
    dfdxs = (A0/(N_T*(xs_E_mon+dxs)*(1-np.exp(-decay_const*t_irr))) - A0/(N_T*(xs_E_mon-dxs)*(1-np.exp(-decay_const*t_irr))))/dxs
    dfdx_list.append(dfdxs)
    unc_list.append(xs_E_mon_unc)

    dt_irr = t_irr*1e-8
    dfdt_irr = (A0/(N_T*xs_E_mon*(1-np.exp(-decay_const*(t_irr+dt_irr)))) - A0/(N_T*xs_E_mon*(1-np.exp(-decay_const*(t_irr-dt_irr)))))/dt_irr 
    dfdx_list.append(dfdt_irr)
    unc_list.append(t_irr_unc)

    dfdx = np.array(dfdx_list)
    unc = np.array(unc_list)
    beam_current_unc = np.sqrt(np.sum(np.multiply(dfdx,dfdx)* np.multiply(unc,unc)))*1.60217634e-19*1e9

    #beam_current_unc = beam_current_in_nA*np.sqrt( (A0_unc/A0)**2 + (N_T_unc/N_T)**2 + (xs_E_mon_unc/xs_E_mon)**2 + (exp_unc/np.exp(-decay_const*t_irr))**2 ) #[nA]

    return beam_current_in_nA, beam_current_unc





def beam_current_list(reaction_list, A0_list, areal_dens_list, Mm_list, decay_const_list, xs_E_mon_list, A0_unc_list, areal_dens_unc_list, Mm_unc_list, xs_E_mon_unc_list):
    beam_current_list = []
    beam_current_unc_list = []

    for i in range(len(reaction_list)):
        beam_current, beam_current_unc = convert_A0_to_beam_current_w_unc(A0_list[i], areal_dens_list[i], Mm_list[i],    
        decay_const_list[i], xs_E_mon_list[i], A0_unc_list[i], areal_dens_unc_list[i], Mm_unc_list[i], xs_E_mon_unc_list[i])
        beam_current_list.append(beam_current)
        beam_current_unc_list.append(beam_current_unc)

    return beam_current_list, beam_current_unc_list