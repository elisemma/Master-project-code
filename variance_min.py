
import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.interpolate import splev, splrep





def calculate_chi_squared(self):
    # Finding the chi squared of the beam current by giving the calculated beam current with uncertainty and the weighted average beam current which will be assigned as the true value
    calculated_beam_currents, beam_currents_unc = self.calculate_beam_current_w_unc()
    weighted_average_beam_current = self.calculate_weighted_average_beam_current()

    weighted_average_beam_current_array = np.zeros(len(calculated_beam_currents))
    weighted_average_beam_current_array.fill(weighted_average_beam_current)
    x_diff = weighted_average_beam_current_array-calculated_beam_currents 

    chi_squared = np.sum( np.multiply(x_diff, x_diff)/ np.multiply(beam_currents_unc, beam_currents_unc))

    return chi_squared



















