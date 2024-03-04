import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 




def caclulate_beam_currents_in_foil(dp, compound):
    stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.3f}.csv')

   

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]

    beam_current_list_of_list = [] #This list will cointain lists of beam currents for all the mon reactions in the foils: [[Ni01:56CO, Ni01:58CO, Ni01:61CU], [Ni02:56CO, Ni02:58CO, Ni02:61CU], ...]
    beam_current_unc_list_of_list = [] #Same as beam_current_list_of_list but with the uncertainties

    beam_energy_in_foil_list_list = []

    reaction_list_list = []

    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']
        # beam_energy_in_foil = row['mu_E']

        A0_by_curie_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_curie.csv')

        if os.path.exists(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv'):
            A0_by_hand_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv')
            A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)

        else:
            A0_concat_df = A0_by_curie_df


        # print('A0_conccat: ', A0_concat_df)

        reaction_list = A0_concat_df['Isotope'].tolist()
        # if compound =='Ni':
        #     reaction_list = ['58CO']
        # if compound == 'Ti':
        #     reaction_list = ['48V']
        # A0_concat_df = A0_concat_df[A0_concat_df['Isotope']==reaction_list[0]]

        A0_list = A0_concat_df['A0'].tolist()
        A0_unc_list = A0_concat_df['A0_unc'].tolist()
        A0_unc_list = [1e18 if np.isinf(value) else value for value in A0_unc_list]
        areal_dens = row['areal_density']

        foil = Foil(foil_name, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, dp)

        foil.assign_molar_mass()
        foil.assign_areal_dens_unc_percent()
        foil.calculate_decay_constant()
        foil.find_monitor_cross_section()
        foil.calculate_beam_currents_w_unc()
        foil.calculate_weighted_average_beam_current()

        foil_beam_cur_list = foil.beam_current_list
        foil_beam_cur_unc_list = foil.beam_current_unc_list
        beam_energy_in_foil = foil.beam_energy_in_foil
        avrg_beam_current = foil.weighted_average_beam_current
        avrg_beam_current_unc = np.sqrt(foil.var_weighted_average_beam_current)

        beam_energy_in_foil_array = np.zeros(len(reaction_list))
        beam_energy_in_foil_array.fill(beam_energy_in_foil)

        beam_current_list_of_list.append(foil_beam_cur_list)
        beam_current_unc_list_of_list.append(foil_beam_cur_unc_list)
        beam_energy_in_foil_list_list.append(beam_energy_in_foil_array)
        reaction_list_list.append(reaction_list)

        print(f'Weighted average beam current for foil {foil_name} is {avrg_beam_current} +- {avrg_beam_current_unc}')

    return beam_current_list_of_list, beam_current_unc_list_of_list, beam_energy_in_foil_list_list, reaction_list_list








#__________________________Running the function________________________________

dp = 0.972
beam_current_list_of_list_Ni, beam_current_unc_list_of_list_Ni, beam_energy_in_foil_list_list_Ni, reaction_list_list_Ni = caclulate_beam_currents_in_foil(dp, 'Ni')
beam_current_list_of_list_Ti, beam_current_unc_list_of_list_Ti, beam_energy_in_foil_list_list_Ti, reaction_list_list_Ti = caclulate_beam_currents_in_foil(dp, 'Ti')




#___________________________Sorting the data___________________________________

# Initialize nested dictionary to store beam currents, uncertainties and energies for each reaction
data_by_reaction = {}

# Iterate through reaction_list_list_Ni
for i, reaction_list in enumerate(reaction_list_list_Ni):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'beam_current': [], 'beam_current_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add beam current, beam_current_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['beam_current'].append(beam_current_list_of_list_Ni[i][j])
        data_by_reaction[reaction]['beam_current_unc'].append(beam_current_unc_list_of_list_Ni[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ni[i][j])


# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'beam_current': [], 'beam_current_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add beam current, beam_current_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['beam_current'].append(beam_current_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['beam_current_unc'].append(beam_current_unc_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti[i][j])



#Scaling the beam current in the Ni05 foil by a factor 12.642
data_by_reaction['58CO']['beam_current'][-1] = data_by_reaction['58CO']['beam_current'][-1]*12.642
data_by_reaction['58CO']['beam_current_unc'][-1] = data_by_reaction['58CO']['beam_current_unc'][-1]*12.642
data_by_reaction['61CU']['beam_current'][-1] = data_by_reaction['61CU']['beam_current'][-1]*12.642
data_by_reaction['61CU']['beam_current_unc'][-1] = data_by_reaction['61CU']['beam_current_unc'][-1]*12.642



# ______________________________Plotting________________________________________

marker_list = ['d', '*', 's', '<', 'o']
color_list = ['deepskyblue', 'mediumseagreen', 'gold', 'violet', 'mediumvioletred']
color_background_list = ['peachpuff', 'lightgreen', 'pink', 'paleturquoise', 'lemonchiffon']

lower_energy_compartments = data_by_reaction['48V']['energy'][:]
upper_energy_compartments = data_by_reaction['58CO']['energy'][:]

lower_energy_compartments.reverse()
upper_energy_compartments.reverse()

compartment_separation_energies = (np.array(upper_energy_compartments[:-1]) + np.array(lower_energy_compartments[1:]))/2



# Loop over every reaction
for i, (reaction, data) in enumerate(data_by_reaction.items()):
    beam_currents = data['beam_current']
    beam_currents_unc = data['beam_current_unc']
    energies = data['energy']

    plt.errorbar(energies, beam_currents, yerr=beam_currents_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=reaction)

# Plot shaded regions
for i in range(len(compartment_separation_energies)+1):
    if i == 0:
        lower_bound = 0
    else:
        lower_bound = compartment_separation_energies[i-1]
    if i == len(compartment_separation_energies):
        upper_bound = upper_energy_compartments[-1]+1
    else:
        upper_bound = compartment_separation_energies[i]

    plt.axvspan(lower_bound, upper_bound, facecolor=color_background_list[i], alpha=0.4)



# plt.xlim([lower_energy_compartments[0]-1, upper_energy_compartments[-1]+1])
plt.xlabel('Beam energy (MeV)')
plt.ylabel('Beam current (nA)')
plt.title(f'dp = {dp:.3f}')
plt.legend()
plt.show()

print(data_by_reaction)


# print('Beam current for 58Co: ', data_by_reaction['58CO'])
# print('Beam current for 61Cu: ', data_by_reaction['61CU'])




