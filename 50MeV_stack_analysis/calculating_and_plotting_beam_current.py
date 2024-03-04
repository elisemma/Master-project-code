import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 

def caclulate_beam_currents_in_foil(dp, compound):
    stack_df = pd.read_csv(f'./Stack_calculations/stack_50MeV_dp_{dp:.3f}.csv')

   

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]
    # monitor_stack_df = monitor_stack_df.drop(monitor_stack_df[monitor_stack_df['name'] == 'Fe01'].index)

    beam_current_list_of_list = [] #This list will cointain lists of beam currents for all the mon reactions in the foils: [[Ti01:46SC, Ti01:48V], [Ti02:46SC, Ti02:48V], ...]
    beam_current_unc_list_of_list = [] #Same as beam_current_list_of_list but with the uncertainties
    beam_energy_in_foil_list_list = []
    reaction_list_list = []

    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']

        A0_by_curie_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_curie.csv')
        A0_concat_df = A0_by_curie_df


        # print('A0_conccat: ', A0_concat_df)

        reaction_list = A0_concat_df['Isotope'].tolist()
        A0_list = A0_concat_df['A0'].tolist()
        A0_unc_list = A0_concat_df['A0_unc'].tolist()
        areal_dens = row['areal_density']

        foil = Foil(foil_name, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, dp)

        foil.assign_molar_mass()
        foil.assign_areal_dens_unc_percent()
        foil.calculate_decay_constant()
        foil.find_monitor_cross_section()
        foil.calculate_beam_currents_w_unc()
        # foil.calculate_weighted_average_beam_current()

        foil_beam_cur_list = foil.beam_current_list
        foil_beam_cur_unc_list = foil.beam_current_unc_list
        beam_energy_in_foil = foil.beam_energy_in_foil
        avrg_beam_current = foil.weighted_average_beam_current

        beam_energy_in_foil_array = np.zeros(len(reaction_list))
        beam_energy_in_foil_array.fill(beam_energy_in_foil)

        beam_current_list_of_list.append(foil_beam_cur_list)
        beam_current_unc_list_of_list.append(foil_beam_cur_unc_list)
        beam_energy_in_foil_list_list.append(beam_energy_in_foil_array)
        reaction_list_list.append(reaction_list)

        # print(f'Weighted average beam current for foil {foil_name} is {avrg_beam_current}')
        # print(foil_beam_cur_list)

    return beam_current_list_of_list, beam_current_unc_list_of_list, beam_energy_in_foil_list_list, reaction_list_list








#______________________Runnig the code____________________________________________
dp = 0.990
beam_current_list_of_list_Fe, beam_current_unc_list_of_list_Fe, beam_energy_in_foil_list_list_Fe, reaction_list_list_Fe = caclulate_beam_currents_in_foil(dp, 'Fe')
beam_current_list_of_list_Ti, beam_current_unc_list_of_list_Ti, beam_energy_in_foil_list_list_Ti, reaction_list_list_Ti = caclulate_beam_currents_in_foil(dp, 'Ti')




#___________________________Sorting the data___________________________________

# Initialize nested dictionary to store beam currents, uncertainties and energies for each reaction
data_by_reaction = {}

# Iterate through reaction_list_list_Ni
for i, reaction_list in enumerate(reaction_list_list_Fe):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'beam_current': [], 'beam_current_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add beam current, beam_current_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['beam_current'].append(beam_current_list_of_list_Fe[i][j])
        data_by_reaction[reaction]['beam_current_unc'].append(beam_current_unc_list_of_list_Fe[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Fe[i][j])


# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'beam_current': [], 'beam_current_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add beam current, beam_current_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['beam_current'].append(beam_current_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['beam_current_unc'].append(beam_current_unc_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti[i][j])


print(data_by_reaction)


# ______________________________Plotting________________________________________

marker_list = ['d', '*', 's']
color_list = ['deepskyblue', 'gold', 'hotpink']
# color_background_list = ['peachpuff', 'lightgreen', 'pink', 'paleturquoise', 'lemonchiffon']

# lower_energy_compartments = data_by_reaction['48V']['energy'][:]
# upper_energy_compartments = data_by_reaction['56CO']['energy'][:]

# lower_energy_compartments.reverse()
# upper_energy_compartments.reverse()

# compartment_separation_energies = (np.array(upper_energy_compartments[:-1]) + np.array(lower_energy_compartments[1:]))/2



# Loop over every reaction
for i, (reaction, data) in enumerate(data_by_reaction.items()):
    beam_currents = data['beam_current']
    beam_currents_unc = data['beam_current_unc']
    energies = data['energy']

    plt.errorbar(energies, beam_currents, yerr=beam_currents_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=reaction)

# # Plot shaded regions
# for i in range(len(compartment_separation_energies)+1):
#     if i == 0:
#         lower_bound = 0
#     else:
#         lower_bound = compartment_separation_energies[i-1]
#     if i == len(compartment_separation_energies):
#         upper_bound = upper_energy_compartments[-1]+1
#     else:
#         upper_bound = compartment_separation_energies[i]

#     plt.axvspan(lower_bound, upper_bound, facecolor=color_background_list[i], alpha=0.4)



# plt.xlim([lower_energy_compartments[0]-1, upper_energy_compartments[-1]+1])
plt.xlabel('Beam energy (MeV)')
plt.ylabel('Beam current (nA)')
plt.title(f'dp = {dp:.3f}')
plt.legend()
plt.show()


















