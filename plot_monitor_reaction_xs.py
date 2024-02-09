import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 




def caclulate_xs_in_foil(dp, compound):
    stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.2f}.csv')
   

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]
    # compartments = ['01', '02', '03', '04']  # Example list of compartments
    # monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}') & stack_df['name'].str.contains('|'.join(compartments))]


    calc_xs_list_of_list = [] #This list will cointain lists of cross sections calculated from beam currents for all the mon reactions in the foils: [[Ni01:56CO, Ni01:58CO, Ni01:61CU], [Ni02:56CO, Ni02:58CO, Ni02:61CU], ...]
    calc_xs_unc_list_of_list = [] #Same as calc_xs_list_of_list but with the uncertainties

    beam_energy_in_foil_list_list = []

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


        # print('A0_conccat: ', A0_concat_df)

        reaction_list = A0_concat_df['Isotope'].tolist()
        A0_list = A0_concat_df['A0'].tolist()
        A0_unc_list = A0_concat_df['A0_unc'].tolist()
        A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]
        areal_dens = row['areal_density']
        areal_dens_unc_percent = 2 #XXXXXXXXXXXX this is not true, need to find it

        foil = Foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent, dp)

        foil.assign_molar_mass()
        foil.calculate_decay_constant()
        foil.find_monitor_cross_section()
        foil.calculate_beam_currents_w_unc()
        foil.convert_beam_current_back_to_xs_w_unc()

        foil_calc_xs_list = foil.calc_xs_list
        foil_calc_xs_unc_list = foil.calc_xs_unc_list


        beam_energy_in_foil_array = np.zeros(len(reaction_list))
        beam_energy_in_foil_array.fill(beam_energy_in_foil)
        beam_energy_in_foil_list_list.append(beam_energy_in_foil_array)
        reaction_list_list.append(reaction_list)
        calc_xs_list_of_list.append(foil_calc_xs_list)
        calc_xs_unc_list_of_list.append(foil_calc_xs_unc_list)
        

    return calc_xs_list_of_list, calc_xs_unc_list_of_list, beam_energy_in_foil_list_list, reaction_list_list






def get_IAEA_monitro_xs(reaction_product):
    filename = './Monitor_cross_section_data/IAEA_monitor_xs_' + reaction_product + '.txt'
    E_mon_list = []
    xs_list = []
    xs_unc_list = []
    with open(filename) as file:
        lines = file.readlines()[7:-1]
        
        for line in lines:
            words = line.split()
            E_mon_list.append(float(words[0]))
            xs_list.append(float(words[1]))
            xs_unc_list.append(float(words[2]))
        file.close()

    return E_mon_list, xs_list, xs_unc_list





#__________________________Running the function________________________________

dp = 1.00
calc_xs_list_of_list_Ni, calc_xs_unc_list_of_list_Ni, beam_energy_in_foil_list_list_Ni, reaction_list_list_Ni = caclulate_xs_in_foil(dp, 'Ni')
calc_xs_list_of_list_Ti, calc_xs_unc_list_of_list_Ti, beam_energy_in_foil_list_list_Ti, reaction_list_list_Ti = caclulate_xs_in_foil(dp, 'Ti')



#___________________________Sorting the data___________________________________

# Initialize nested dictionary to store beam currents, uncertainties and energies for each reaction
data_by_reaction = {}

# Iterate through reaction_list_list_Ni
for i, reaction_list in enumerate(reaction_list_list_Ni):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list_Ni[i][j])
        data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ni[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ni[i][j])



# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti[i][j])





#____________________________Getting monitor xs__________________________________


E_mon_list_46SC, xs_list_46SC, xs_unc_list_46SC = get_IAEA_monitro_xs('46SC')
E_mon_list_48V, xs_list_48V, xs_unc_list_48V = get_IAEA_monitro_xs('48V')
E_mon_list_56CO, xs_list_56CO, xs_unc_list_56CO = get_IAEA_monitro_xs('56CO')
E_mon_list_58CO, xs_list_58CO, xs_unc_list_58CO = get_IAEA_monitro_xs('58CO')
E_mon_list_61CU, xs_list_61CU, xs_unc_list_61CU = get_IAEA_monitro_xs('61CU')






# ______________________________Plotting________________________________________

marker_list = ['d', '*', 's', '<', 'o']
color_list = ['deepskyblue', 'mediumseagreen', 'gold', 'violet', 'mediumvioletred']


# Loop over every reaction
for i, (reaction, data) in enumerate(data_by_reaction.items()):
    calc_xs = data['calc_xs']
    calc_xs_unc = data['calc_xs_unc']
    energies = data['energy']

    plt.errorbar(energies, calc_xs, yerr=calc_xs_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=reaction)
plt.plot(E_mon_list_56CO, xs_list_56CO, color='deepskyblue', label='IAEA 56Co')
plt.plot(E_mon_list_58CO, xs_list_58CO, color='mediumseagreen', label='IAEA 58Co')
plt.plot(E_mon_list_61CU, xs_list_61CU, color='gold', label='IAEA 61Cu')
plt.plot(E_mon_list_48V, xs_list_48V, color='violet', label='IAEA 48V')
plt.plot(E_mon_list_46SC, xs_list_46SC, color='mediumvioletred', label='IAEA 46Sc')

plt.xlabel('Beam energy (MeV)')
plt.ylabel('Cross section (mb)')
plt.title(f'dp = {dp:.2f}, beam current = 126nA for all foils')
plt.legend()
plt.show()

# plt.plot(E_mon_list_56CO, xs_list_56CO, color='deepskyblue', label='IAEA 56Co')
# plt.plot(E_mon_list_58CO, xs_list_58CO, color='mediumseagreen', label='IAEA 58Co')
# plt.plot(E_mon_list_61CU, xs_list_61CU, color='gold', label='IAEA 61Cu')
# plt.plot(E_mon_list_48V, xs_list_48V, color='violet', label='IAEA 48V')
# plt.plot(E_mon_list_46SC, xs_list_46SC, color='mediumvioletred', label='IAEA 46Sc')

# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'dp = {dp:.2f}')
# plt.legend()
# plt.show()






