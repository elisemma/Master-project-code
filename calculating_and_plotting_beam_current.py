import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil



beam_energy_in_foil = 27.544845
target_material = 'Ni'
reaction_list = ['56CO', '58CO', '61CU']
A0_list = [256.99, 2850, 220000]
A0_unc_list = [20, 300, 20000]
areal_dens = 23.253 #[g/cm^2]
areal_dens_unc_percent = 2

Ni01 = Foil(beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent)

Ni01.assign_molar_mass()
Ni01.calculate_decay_constant()
Ni01.find_monitor_cross_section()
Ni01.calculate_beam_currents_w_unc()
Ni01.calculate_weighted_average_beam_current()


Ni01_average_beam_cur = Ni01.weighted_average_beam_current
Ni01_average_beam_cur_var = Ni01.var_weighted_average_beam_current

print(f'Avr beam current in foil Ni01 = {Ni01_average_beam_cur} +- {np.sqrt(Ni01_average_beam_cur_var)}')


Ni01_beam_cur_list = Ni01.beam_current_list
Ni01_beam_cur_unc_list = Ni01.beam_current_unc_list

for i in range(len(Ni01_beam_cur_list)):
    plt.errorbar(beam_energy_in_foil, Ni01_beam_cur_list[i], yerr = Ni01_beam_cur_unc_list[i], label = reaction_list[i])
plt.errorbar(beam_energy_in_foil, Ni01_average_beam_cur, yerr = np.sqrt(Ni01_average_beam_cur_var), label = 'Average')
plt.legend()
plt.show()
