import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 



def caclulate_xs_in_foil(dp, compound):
    stack_df = pd.read_csv(f'./Stack_calculations/stack_50MeV_dp_{dp:.2f}.csv')
    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]
    monitor_stack_df = monitor_stack_df.drop(monitor_stack_df[monitor_stack_df['name'] == 'Fe01'].index)


    calc_xs_list_of_list = [] 
    calc_xs_unc_list_of_list = []

    beam_energy_in_foil_list_list = []

    reaction_list_list = []


    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']
        # beam_energy_in_foil = row['mu_E']

        A0_by_curie_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_curie.csv')
        A0_concat_df = A0_by_curie_df

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
        foil.calculate_weighted_average_beam_current()
        foil.convert_beam_current_back_to_xs_w_unc()

        foil_calc_xs_list = foil.calc_xs_list
        foil_calc_xs_unc_list = foil.calc_xs_unc_list
        beam_energy_in_foil = foil.beam_energy_in_foil


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



# #_________________________Energy info from another file________________________
# foil_energy_data = {'Ni01': {'energy': 27.337003000000003, 'min_unc': 0.6070030000000024, 'plus_unc': 0.6529969999999956}, 
#                     'Zr01': {'energy': 26.375581480000005, 'min_unc': 0.605581480000005, 'plus_unc': 0.5944185199999943}, 
#                     'Ti01': {'energy': 25.518248420000003, 'min_unc': 0.5882484200000029, 'plus_unc': 0.6117515799999964}, 
#                     'Ni02': {'energy': 20.392268865, 'min_unc': 0.7422688649999998, 'plus_unc': 0.7577311350000002}, 
#                     'Zr02': {'energy': 19.1745053, 'min_unc': 0.7845052999999993, 'plus_unc': 0.7754946999999994}, 
#                     'Ti02': {'energy': 18.058134579999997, 'min_unc': 0.8081345799999973, 'plus_unc': 0.8118654200000037}, 
#                     'Ni03': {'energy': 13.9680305, 'min_unc': 0.9780305000000009, 'plus_unc': 1.0619695}, 
#                     'Zr03': {'energy': 12.309193775999999, 'min_unc': 1.059193775999999, 'plus_unc': 1.1008062240000012}, 
#                     'Ti03': {'energy': 10.706294196000002, 'min_unc': 1.1362941960000015, 'plus_unc': 1.2037058039999984}, 
#                     'Ni04': {'energy': 8.648899379944393, 'min_unc': 1.2988993799443929, 'plus_unc': 1.5211006200556056}, 
#                     'Zr04': {'energy': 6.082350534817555, 'min_unc': 1.4923505348175548, 'plus_unc': 1.8676494651824447}, 
#                     'Ti04': {'energy': 3.6998866038357203, 'min_unc': 2.2298866038357206, 'plus_unc': 2.090113396164279}, 
#                     'Ni05': {'energy': 2.2266116583208797, 'min_unc': 1.8966116583208799, 'plus_unc': 1.4033883416791202}, 
#                     'Zr05': {'energy': 1.5675066930508672, 'min_unc': 1.2375066930508674, 'plus_unc': 0.7424933069491328}, 
#                     'Ti05': {'energy': 0.7071428571428571, 'min_unc': 0.37714285714285717, 'plus_unc': 0.2828571428571429}}




#__________________________Running the function________________________________

dp = 0.97
calc_xs_list_of_list_Fe, calc_xs_unc_list_of_list_Fe, beam_energy_in_foil_list_list_Fe, reaction_list_list_Fe = caclulate_xs_in_foil(dp, 'Fe')
calc_xs_list_of_list_Ti, calc_xs_unc_list_of_list_Ti, beam_energy_in_foil_list_list_Ti, reaction_list_list_Ti = caclulate_xs_in_foil(dp, 'Ti')


# calc_xs_list_of_list_Fe_p0, calc_xs_unc_list_of_list_Fe_p0, beam_energy_in_foil_list_list_Fe_p0, reaction_list_list_Fe_p0 = caclulate_xs_in_foil(0.988, 'Fe')
# calc_xs_list_of_list_Ti_p0, calc_xs_unc_list_of_list_Ti_p0, beam_energy_in_foil_list_list_Ti_p0, reaction_list_list_Ti_p0 = caclulate_xs_in_foil(0.988, 'Ti')

# calc_xs_list_of_list_Fe_p1, calc_xs_unc_list_of_list_Fe_p1, beam_energy_in_foil_list_list_Fe_p1, reaction_list_list_Fe_p1 = caclulate_xs_in_foil(0.977, 'Fe')
# calc_xs_list_of_list_Ti_p1, calc_xs_unc_list_of_list_Ti_p1, beam_energy_in_foil_list_list_Ti_p1, reaction_list_list_Ti_p1 = caclulate_xs_in_foil(0.977, 'Ti')




#___________________________Sorting the data___________________________________

# Initialize nested dictionary to store beam currents, uncertainties and energies for each reaction
data_by_reaction = {}

# Iterate through reaction_list_list_Fe
for i, reaction_list in enumerate(reaction_list_list_Fe):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list_Fe[i][j])
        data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Fe[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Fe[i][j])



# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction:
            data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti[i][j])
        data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti[i][j])


# data_by_reaction_p0 = {}

# # Iterate through reaction_list_list_Ni
# for i, reaction_list in enumerate(reaction_list_list_Ni_p0):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_p0:
#             data_by_reaction_p0[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_p0[reaction]['calc_xs'].append(calc_xs_list_of_list_Ni_p0[i][j])
#         data_by_reaction_p0[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ni_p0[i][j])
#         data_by_reaction_p0[reaction]['energy'].append(beam_energy_in_foil_list_list_Ni_p0[i][j])



# # Iterate through reaction_list_list_Ti
# for i, reaction_list in enumerate(reaction_list_list_Ti_p0):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_p0:
#             data_by_reaction_p0[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_p0[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_p0[i][j])
#         data_by_reaction_p0[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_p0[i][j])
#         data_by_reaction_p0[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_p0[i][j])


# data_by_reaction_p1 = {}

# # Iterate through reaction_list_list_Ni
# for i, reaction_list in enumerate(reaction_list_list_Ni_p1):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_p1:
#             data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Ni_p1[i][j])
#         data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ni_p1[i][j])
#         data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Ni_p1[i][j])



# # Iterate through reaction_list_list_Ti
# for i, reaction_list in enumerate(reaction_list_list_Ti_p1):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_p1:
#             data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_p1[i][j])
#         data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_p1[i][j])
#         data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_p1[i][j])






#____________________________Getting monitor xs__________________________________

E_mon_list_46SC, xs_list_46SC, xs_unc_list_46SC = get_IAEA_monitro_xs('46SC')
E_mon_list_48V, xs_list_48V, xs_unc_list_48V = get_IAEA_monitro_xs('48V')
E_mon_list_56CO, xs_list_56CO, xs_unc_list_56CO = get_IAEA_monitro_xs('56CO')






# ______________________________Plotting________________________________________

marker_list = ['d', '*', 's']
color_list = ['deepskyblue','gold', 'hotpink']


# Loop over every reaction
for i, (reaction, data) in enumerate(data_by_reaction.items()):
    calc_xs = data['calc_xs']
    calc_xs_unc = data['calc_xs_unc']
    energies = data['energy']
    plt.errorbar(energies, calc_xs, yerr=calc_xs_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=f'{reaction} calculated with avrg bc')

plt.plot(E_mon_list_56CO, xs_list_56CO, color='deepskyblue', label='IAEA 56Co')
plt.plot(E_mon_list_46SC, xs_list_46SC, color='gold', label='IAEA 46Sc')
plt.plot(E_mon_list_48V, xs_list_48V, color='hotpink', label='IAEA 48V')


plt.xlabel('Beam energy (MeV)')
plt.ylabel('Cross section (mb)')
plt.title(f'dp = {dp:.2f}')
plt.legend()
plt.show()

# Fe_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Fe')]
# Fe_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Fe')]
# Ti_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Ti')]
# Ti_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Ti')]


# calc_xs_56CO_p0 = data_by_reaction_p0['56CO']['calc_xs']
# calc_xs_unc_56CO_p0 = data_by_reaction_p0['56CO']['calc_xs_unc']
# energies_56CO_p0 = data_by_reaction_p0['56CO']['energy']

# calc_xs_46SC_p0 = data_by_reaction_p0['46SC']['calc_xs']
# calc_xs_unc_46SC_p0 = data_by_reaction_p0['46SC']['calc_xs_unc']
# energies_46SC_p0 = data_by_reaction_p0['46SC']['energy']

# calc_xs_48V_p0 = data_by_reaction_p0['48V']['calc_xs']
# calc_xs_unc_48V_p0 = data_by_reaction_p0['48V']['calc_xs_unc']
# energies_48V_p0 = data_by_reaction_p0['48V']['energy']


# calc_xs_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs']
# calc_xs_unc_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs_unc']
# energies_56CO_p1 = data_by_reaction_p1['56CO']['energy']

# calc_xs_46SC_p1 = data_by_reaction_p1['46SC']['calc_xs']
# calc_xs_unc_46SC_p1 = data_by_reaction_p1['46SC']['calc_xs_unc']
# energies_46SC_p1 = data_by_reaction_p1['46SC']['energy']

# calc_xs_48V_p1 = data_by_reaction_p1['48V']['calc_xs']
# calc_xs_unc_48V_p1 = data_by_reaction_p1['48V']['calc_xs_unc']
# energies_48V_p1 = data_by_reaction_p1['48V']['energy']


# plt.errorbar(energies_56CO_p0, calc_xs_56CO_p0, xerr=[Fe_energy_min_unc_list[:-1], Fe_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
# plt.errorbar(energies_56CO_p1, calc_xs_56CO_p1, xerr=[Fe_energy_min_unc_list[:-1], Fe_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
# plt.plot(E_mon_list_56CO, xs_list_56CO, color='hotpink', label='IAEA')
# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'56CO')
# plt.legend()
# plt.xlim(0,50)
# plt.show()

# plt.errorbar(energies_46SC_p0, calc_xs_46SC_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
# plt.errorbar(energies_46SC_p1, calc_xs_46SC_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
# plt.plot(E_mon_list_46SC, xs_list_46SC, color='hotpink', label='IAEA')
# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'46SC')
# plt.legend()
# plt.xlim(0,50)
# plt.show()

# plt.errorbar(energies_48V_p0, calc_xs_48V_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
# plt.errorbar(energies_48V_p1, calc_xs_48V_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
# plt.plot(E_mon_list_48V, xs_list_48V, color='hotpink', label='IAEA')
# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'48V')
# plt.legend()
# plt.xlim(0,50)
# plt.show()





