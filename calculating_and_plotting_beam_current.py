import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 




def caclulate_beam_currents_in_foil(dp, compound):
    # stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.2f}.csv')
    stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations_copy_dp_080_120/stack_30MeV_dp_{dp:.2f}.csv')
    

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]
    # compartment = '03'
    # monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')  & (stack_df['name'].str.contains(compartment))]

    beam_current_list_of_list = [] #This list will cointain lists of beam currents for all the mon reactions in the foils: [[Ni01:56CO, Ni01:58CO, Ni01:61CU], [Ni02:56CO, Ni02:58CO, Ni02:61CU], ...]
    beam_current_unc_list_of_list = [] #Same as beam_current_list_of_list but with the uncertainties

    beam_energy_in_foil_list_list = []

    reaction_list_list = []

    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']
        beam_energy_in_foil = row['mu_E']

        A0_by_curie_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_curie.csv')

        if os.path.exists(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv'):
            A0_by_hand_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv')
            A0_by_hand_second_ord_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand_second_order.csv')


            A0_concat_df = pd.concat((A0_by_hand_second_ord_df, A0_by_hand_df, A0_by_curie_df), axis=0)
        else:
            A0_concat_df = A0_by_curie_df


        # print('A0_conccat: ', A0_concat_df)

        reaction_list = A0_concat_df['Isotope'].tolist()
        A0_list = A0_concat_df['A0'].tolist()
        A0_unc_list = A0_concat_df['A0_unc'].tolist()
        A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]
        areal_dens = row['areal_density']
        areal_dens_unc_percent = 2 #XXXXXXXXXXXX this is not true, need to find it

        foil = Foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent, dp=1.00)

        foil.assign_molar_mass()
        foil.calculate_decay_constant()
        foil.find_monitor_cross_section()
        foil.calculate_beam_currents_w_unc()
        foil.calculate_weighted_average_beam_current()


        foil_average_beam_cur = foil.weighted_average_beam_current
        foil_average_beam_cur_var = foil.var_weighted_average_beam_current
        foil_beam_cur_list = foil.beam_current_list
        foil_beam_cur_unc_list = foil.beam_current_unc_list

        beam_energy_in_foil_array = np.zeros(len(reaction_list))
        beam_energy_in_foil_array.fill(beam_energy_in_foil)

        beam_current_list_of_list.append(foil_beam_cur_list)
        beam_current_unc_list_of_list.append(foil_beam_cur_unc_list)
        beam_energy_in_foil_list_list.append(beam_energy_in_foil_array)
        reaction_list_list.append(reaction_list)


    return beam_current_list_of_list, beam_current_unc_list_of_list, beam_energy_in_foil_list_list, reaction_list_list










#__________________________Running the function________________________________

dp = 1.00
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




# ______________________________Plotting________________________________________

marker_list = ['d', '*', 's', '<', 'o']
color_list = ['deepskyblue', 'mediumseagreen', 'gold', 'violet', 'mediumvioletred']
color_background_list = ['peachpuff', 'lightgreen', 'pink', 'paleturquoise', 'lemonchiffon']

lower_energy_compartments = data_by_reaction['48V']['energy']
upper_energy_compartments = data_by_reaction['61CU']['energy']

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

    plt.axvspan(lower_bound, upper_bound, facecolor=color_background_list[i], alpha=0.4)#, label=f'compartment {5-i}')


plt.xlim([lower_energy_compartments[0]-1, upper_energy_compartments[-1]+1])
plt.xlabel('Beam energy (MeV)')
plt.ylabel('Beam current (nA)')
plt.title(f'dp = {dp:.2f}')
plt.legend()
plt.show()






