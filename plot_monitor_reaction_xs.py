import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 



print('comp 5')

def caclulate_xs_in_foil(dp, compound):
    stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.2f}.csv')
   

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]
    compartments = ['05']  # Example list of compartments
    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}') & stack_df['name'].str.contains('|'.join(compartments))]


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
        lines = file.readlines()[6:-1]
        
        for line in lines:
            words = line.split()
            E_mon_list.append(float(words[0]))
            xs_list.append(float(words[1]))
            xs_unc_list.append(float(words[2]))
        file.close()

    return E_mon_list, xs_list, xs_unc_list





#__________________________Running the function________________________________

dp = 0.98
# calc_xs_list_of_list_Ni, calc_xs_unc_list_of_list_Ni, beam_energy_in_foil_list_list_Ni, reaction_list_list_Ni = caclulate_xs_in_foil(dp, 'Ni')
# calc_xs_list_of_list_Ti, calc_xs_unc_list_of_list_Ti, beam_energy_in_foil_list_list_Ti, reaction_list_list_Ti = caclulate_xs_in_foil(dp, 'Ti')



# #___________________________Sorting the data___________________________________

# # Initialize nested dictionary to store beam currents, uncertainties and energies for each reaction
# data_by_reaction = {}

# # Iterate through reaction_list_list_Ni
# for i, reaction_list in enumerate(reaction_list_list_Ni):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction:
#             data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list_Ni[i][j])
#         data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ni[i][j])
#         data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ni[i][j])



# # Iterate through reaction_list_list_Ti
# for i, reaction_list in enumerate(reaction_list_list_Ti):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction:
#             data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti[i][j])
#         data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti[i][j])
#         data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti[i][j])



# data_by_reaction = { #dp = 1.00
#     '56CO': {
#         'calc_xs': [10.094852830778155, 8.15643300328985, 20.977958528132273, 24.688210900077042],
#         'calc_xs_unc': [4.594342020913995e-37, 4.8524957987970285e-37, 3.9109370159315166e-37, 3.755572816188005e-37],
#         'energy': [27.27237216, 20.12526522, 13.427759768571429, 7.682573316]
#     },
#     '58CO': {
#         'calc_xs': [132.16792298526033, 218.63067586236684, 177.0206591262901, 38.582133229951005, 0.8845961247537194],
#         'calc_xs_unc': [5.6676530016896266e-33, 6.1351635502271986e-33, 4.705708136783849e-37, 5.6452377974747526e-33, 5.568309880475e-33],
#         'energy': [27.27237216, 20.12526522, 13.427759768571429, 7.682573316, 0.0798290376470588]
#     },
#     '61CU': {
#         'calc_xs': [15.840902686067421, 16.399755118277493, 26.98457947266417, 50.1854945505117, 0.03261116248742061],
#         'calc_xs_unc': [9.278557392630717e-40, 6.632468418604548e-40, 7.064580257226981e-40, 4.271571867077405e-40, 1.1314528671316037e-35],
#         'energy': [27.27237216, 20.12526522, 13.427759768571429, 7.682573316, 0.0798290376470588]
#     },
#     '48V': {
#         'calc_xs': [186.85609649561695, 330.14823990162216, 121.71669096533111, 3.5178926785055236],
#         'calc_xs_unc': [1.5014215900481878e-37, 1.5049285299293398e-37, 1.848158415386213e-37, 2.364270643308602e-37],
#         'energy': [25.40594876, 17.70611474, 9.9478172, 1.336401931764706]
#     },
#     '46SC': {
#         'calc_xs': [23.8833945989943, 25.590800592954498, 26.78783421726758, 0.632178384575723],
#         'calc_xs_unc': [7.518923231909391e-37, 7.236804758769653e-37, 6.24454676546448e-37, 4.556866750696946e-37],
#         'energy': [25.40594876, 17.70611474, 9.9478172, 1.336401931764706]
#     }
# }



data_by_reaction = { #dp=0.98
    '56CO': {
        'calc_xs': [10.10437739072342, 8.060189980645456, 20.53545636596296, 29.02932121790106],
        'calc_xs_unc': [4.533286962237811e-37, 4.928800131787287e-37, 4.42586156218653e-37, 4.12746418842812e-37],
        'energy': [27.32960846, 20.358895125, 13.90021321, 8.530354149230767]
    },
    '58CO': {
        'calc_xs': [132.29262428862913, 216.0509112668269, 173.28664352572662, 45.36631444608274, 1.294133827783704],
        'calc_xs_unc': [5.673000469658724e-33, 6.062770791765247e-33, 3.911165600709885e-37, 6.6378815889761e-33, 8.146232317113121e-33],
        'energy': [27.32960846, 20.358895125, 13.90021321, 8.530354149230767, 0.2934494952941176]
    },
    '61CU': {
        'calc_xs': [15.855848681789169, 16.206243812223207, 26.415375622543557, 59.00998041868967, 0.0477090135908879],
        'calc_xs_unc': [9.320458478936302e-40, 6.593881191617527e-40, 7.591847222398609e-40, 6.356964092117998e-40, 1.6552767855550727e-35],
        'energy': [27.32960846, 20.358895125, 13.90021321, 8.530354149230767, 0.2934494952941176]
    },
    '48V': {
        'calc_xs': [186.04241895493107, 326.2372098921215, 145.38603195821898, 8.3753688520458, 0.003397971540980206],
        'calc_xs_unc': [1.516922013517768e-37, 1.4899419357185178e-37, 1.5433561095606713e-37, 2.6349808663657038e-37, 1.5852341419673778e-37],
        'energy': [25.50470122, 18.013976199999995, 10.611885444, 2.674746648947369, 0.0300008647058823]
    },
    '46SC': {
        'calc_xs': [23.779392737963253, 25.287644686032145, 31.997065404221676, 1.5050849002481383, 0.0030083104101382655],
        'calc_xs_unc': [7.439882323308149e-37, 7.31246688758761e-37, 7.657673542612117e-37, 3.7578444666618648e-37, 2.7245753048984864e-34],
        'energy': [25.50470122, 18.013976199999995, 10.611885444, 2.674746648947369, 0.0300008647058823]
    }
}

print(data_by_reaction)






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
    plt.errorbar(energies, calc_xs, yerr=calc_xs_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=f'{reaction} calculated with avrg bc')

plt.plot(E_mon_list_56CO, xs_list_56CO, color='deepskyblue', label='IAEA 56Co')
plt.plot(E_mon_list_58CO, xs_list_58CO, color='mediumseagreen', label='IAEA 58Co')
plt.plot(E_mon_list_61CU, xs_list_61CU, color='gold', label='IAEA 61Cu')
plt.plot(E_mon_list_48V, xs_list_48V, color='violet', label='IAEA 48V')
plt.plot(E_mon_list_46SC, xs_list_46SC, color='mediumvioletred', label='IAEA 46Sc')

plt.xlabel('Beam energy (MeV)')
plt.ylabel('Cross section (mb)')
plt.title(f'dp = {dp:.2f}')
plt.legend()
plt.show()









