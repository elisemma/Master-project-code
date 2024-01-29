
import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.optimize import curve_fit


# Plan: Get the beam currents with uncertainty for one compartment
# Then make a p0 fit whith these data points
# Give the beam currents as observed and the p0 fit as the true value to chi2 
# plot the chi2 as a function of dp and find the minimum 

def p0(x, a):
    return a


def fit_p0(data, unc_data):
    x_array = np.arange(len(data))
    popt, cov  = curve_fit(p0, x_array, data, p0=100, sigma=unc_data)
    return popt[0]


def calculate_chi2(observed, expected):
    diff = observed-expected 
    chi2 = np.sum(np.multiply(diff, diff)/expected)
    return chi2


def run_chi2(data, unc_data):
    true = fit_p0(data, unc_data)
    true_array = np.zeros(len(data))
    true_array.fill(true)
    chi2 = calculate_chi2(data, true_array)

    print(chi2)

    plt.errorbar(np.arange(len(data)), data, yerr=unc_data)
    plt.plot(np.arange(len(data)), true_array)
    plt.text(0.05, 0.95, f'chi2 = {chi2}', transform=plt.gca().transAxes, fontsize=12,verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    plt.show()
    return chi2



def plot_chi2(dp_list):
    chi2_list = []
    for dp in dp_list:
        # calculate beam current here
        chi2 = run_chi2(beam_current_list,beam_current_unc_list)
        chi2_list.append(chi2)


    plt.plot(dp, chi2)
    plt.show()
    









test_data = np.array([100, 120, 110, 117, 122])
test_data_unc = np.array([8, 12, 10, 17, 12])


run_the_code(test_data, test_data_unc)














# def calculate_chi_squared():
#     # Finding the chi squared of the beam current by giving the calculated beam current with uncertainty and the weighted average beam current which will be assigned as the true value
#     calculated_beam_currents, beam_currents_unc = self.calculate_beam_current_w_unc()
#     weighted_average_beam_current = self.calculate_weighted_average_beam_current()

#     weighted_average_beam_current_array = np.zeros(len(calculated_beam_currents))
#     weighted_average_beam_current_array.fill(weighted_average_beam_current)
#     x_diff = weighted_average_beam_current_array-calculated_beam_currents 

#     # chi_squared = np.sum( np.multiply(x_diff, x_diff)/ np.multiply(beam_currents_unc, beam_currents_unc))

#     chi_squared = np.sum(np.multiply(x_diff, x_diff)/true)

#     return chi_squared



















