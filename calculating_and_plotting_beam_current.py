import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 





stack_df = pd.read_csv('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_1.00.csv')

monitor_stack_df = stack_df[stack_df['name'].str.contains('Ni|Ti')]
print(monitor_stack_df)

for index, row in monitor_stack_df.iterrows():
    foil_name = row['name']
    # if 'Ni' in foil_name:
    #     reaction_list = ['56CO', '61CU']
    # if 'Ti' in foil_name:
    #     reaction_list = ['46SC', '48V']
    target_material = row['compound']
    beam_energy_in_foil = row['mu_E']

    A0_by_curie_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_curie.csv')

    if os.path.exists(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv'):
        A0_by_hand_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv')

        A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)
    else:
        A0_concat_df = A0_by_curie_df

    print(A0_concat_df)

    reaction_list = A0_concat_df['Isotope'].tolist()
    A0_list = A0_concat_df['A0'].tolist()
    A0_unc_list = A0_concat_df['A0_unc'].tolist()
    A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]
    areal_dens = row['areal_density']
    areal_dens_unc_percent = 2 #XXXXXXXXXXXX this is not true, need to find it
    

    print('reaction_list', reaction_list)
    print('A0', A0_list)
    print('A0_unc', A0_unc_list)

    foil = Foil(beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent)

    foil.assign_molar_mass()
    foil.calculate_decay_constant()
    foil.find_monitor_cross_section()
    foil.calculate_beam_currents_w_unc()
    foil.calculate_weighted_average_beam_current()


    foil_average_beam_cur = foil.weighted_average_beam_current
    foil_average_beam_cur_var = foil.var_weighted_average_beam_current
    foil_beam_cur_list = foil.beam_current_list
    foil_beam_cur_unc_list = foil.beam_current_unc_list

    for i in range(len(foil_beam_cur_list)):
        plt.errorbar(beam_energy_in_foil, foil_beam_cur_list[i], yerr = foil_beam_cur_unc_list[i], label = f'{foil_name}: {reaction_list[i]}')
    plt.errorbar(beam_energy_in_foil, foil_average_beam_cur, yerr = np.sqrt(foil_average_beam_cur_var), label = f'{foil_name}: Average')
# plt.legend()
plt.show()

