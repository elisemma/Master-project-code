
import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.interpolate import splev, splrep





class VarianceMin:
    def __init__(self, stack, beam_energy, projectile, target_material, reaction_list, A0, areal_dens, areal_dens_unc_percent):
        self.stack = stack
        self.beam_energy = beam_energy
        self.projectile = projectile
        self.reaction_list = reaction_list
        self.target_material = target_material
        self.A0 = A0
        self.areal_dens = areal_dens
        self.areal_dens_unc_percent = areal_dens_unc_percent
        self.molar_mass = None
        self.molar_mass_unc = None
        self.decay_const_list = [] 
        

    def energy_in_foil(self, N, E0, dp): 
        # Using the curie data package and the stack information to calculate the energy in all of the foils in the stack
        st = ci.Stack(self.stack, E0=E0, N=N, particle=self.projectile, dp = dp)
        st.saveas(f'stack_calc_E_{E0}_dp_{dp}.csv')
        st_pd = pd.read_csv(f'stack_calc_E_{E0}_dp_{dp}.csv')

        return st_pd




    def calculate_chi_squared(self):
        # Finding the chi squared of the beam current by giving the calculated beam current with uncertainty and the weighted average beam current which will be assigned as the true value
        calculated_beam_currents, beam_currents_unc = self.calculate_beam_current_w_unc()
        weighted_average_beam_current = self.calculate_weighted_average_beam_current()

        weighted_average_beam_current_array = np.zeros(len(calculated_beam_currents))
        weighted_average_beam_current_array.fill(weighted_average_beam_current)
        x_diff = weighted_average_beam_current_array-calculated_beam_currents 

        chi_squared = np.sum( np.multiply(x_diff, x_diff)/ np.multiply(beam_currents_unc, beam_currents_unc))

        return chi_squared


    def find_optimal_beam_current(self):
        # Varying the beam energy and areal densities of the targets in the stack to get the optimal beam current by using the chi squared minimization method























