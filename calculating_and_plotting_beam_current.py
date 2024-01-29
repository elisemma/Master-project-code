import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 




def caclulate_beam_currents_in_foil(dp, compound):
    stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.2f}.csv')

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]

    beam_current_list_of_list = [] #This list will cointain lists of beam currents for all the mon reactions in the foils: [[Ni01:56CO, Ni01:58CO, Ni01:61CU], [Ni02:56CO, Ni02:58CO, Ni02:61CU], ...]
    beam_current_unc_list_of_list = [] #Same as beam_current_list_of_list but with the uncertainties

    beam_energy_in_foil_list = []

    reaction_list_list = []

    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']
        beam_energy_in_foil = row['mu_E']

        A0_by_curie_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_curie.csv')

        if os.path.exists(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv'):
            A0_by_hand_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv')

            A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)
        else:
            A0_concat_df = A0_by_curie_df

        reaction_list = A0_concat_df['Isotope'].tolist()
        A0_list = A0_concat_df['A0'].tolist()
        A0_unc_list = A0_concat_df['A0_unc'].tolist()
        A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]
        areal_dens = row['areal_density']
        areal_dens_unc_percent = 2 #XXXXXXXXXXXX this is not true, need to find it

        foil = Foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent)

        foil.assign_molar_mass()
        foil.calculate_decay_constant()
        foil.find_monitor_cross_section()
        foil.calculate_beam_currents_w_unc()
        foil.calculate_weighted_average_beam_current()


        foil_average_beam_cur = foil.weighted_average_beam_current
        foil_average_beam_cur_var = foil.var_weighted_average_beam_current
        foil_beam_cur_list = foil.beam_current_list
        foil_beam_cur_unc_list = foil.beam_current_unc_list

        beam_current_list_of_list.append(foil_beam_cur_list)
        beam_current_unc_list_of_list.append(foil_beam_cur_unc_list)
        beam_energy_in_foil_list.append(beam_energy_in_foil)
        reaction_list_list.append(reaction_list)


    return beam_current_list_of_list, beam_current_unc_list_of_list, beam_energy_in_foil_list, reaction_list_list



#GOAL: get one list of beam currents for each reaction


beam_current_list_of_list_Ni, beam_current_unc_list_of_list_Ni, beam_energy_in_foil_list_Ni, reaction_list_list_Ni = caclulate_beam_currents_in_foil(1.00, 'Ni')
beam_current_list_of_list_Ti, beam_current_unc_list_of_list_Ti, beam_energy_in_foil_list_Ti, reaction_list_list_Ti = caclulate_beam_currents_in_foil(1.00, 'Ti')



marker_list_Ni = ['d', '*', 's']
marker_list_Ti = ['<', 'o']
color_list_Ni = ['royalblue', 'mediumseagreen', 'gold']
color_list_Ti = ['violet', 'mediumvioletred']

for beam_current_list, beam_current_unc_list, beam_energy, reaction_list in zip(beam_current_list_of_list_Ni, beam_current_unc_list_of_list_Ni, beam_energy_in_foil_list_Ni, reaction_list_list_Ni):

    for i, beam_current, beam_current_unc, reaction in zip(range(len(reaction_list)), beam_current_list, beam_current_unc_list, reaction_list):

        plt.errorbar(beam_energy, beam_current, yerr=beam_current_unc, marker=marker_list_Ni[i], markersize=5, color=color_list_Ni[i], label=reaction)
    

for beam_current_list, beam_current_unc_list, beam_energy, reaction_list in zip(beam_current_list_of_list_Ti, beam_current_unc_list_of_list_Ti, beam_energy_in_foil_list_Ti, reaction_list_list_Ti):

    for i, beam_current, beam_current_unc, reaction in zip(range(len(reaction_list)), beam_current_list, beam_current_unc_list, reaction_list):

        plt.errorbar(beam_energy, beam_current, yerr=beam_current_unc, marker=marker_list_Ti[i], markersize=5, color=color_list_Ti[i], label=reaction)

plt.legend()
plt.show()




