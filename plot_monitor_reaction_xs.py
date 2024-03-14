import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 
from scipy.interpolate import PchipInterpolator




def caclulate_xs_in_foil(dp, compound):
    stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.3f}.csv')

   

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]

    calc_xs_list_of_list = [] #This list will cointain lists of cross sections calculated from beam currents for all the mon reactions in the foils: [[Ni01:56CO, Ni01:58CO, Ni01:61CU], [Ni02:56CO, Ni02:58CO, Ni02:61CU], ...]
    calc_xs_unc_list_of_list = [] #Same as calc_xs_list_of_list but with the uncertainties

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
        A0_list = A0_concat_df['A0'].tolist()
        A0_unc_list = A0_concat_df['A0_unc'].tolist()
        A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]
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
        weighted_average_beam_current = foil.weighted_average_beam_current
        var_weighted_average_beam_current = float(foil.var_weighted_average_beam_current)
        sig_avrg_bc = np.sqrt(var_weighted_average_beam_current)

        print(f'Foil: {foil_name}, avrg bc = {weighted_average_beam_current} +- {sig_avrg_bc}')


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
    interp_xs = PchipInterpolator(E_mon_list, xs_list)
    interp_unc_xs = PchipInterpolator(E_mon_list, xs_unc_list)

    # return E_mon_list, xs_list, xs_unc_list
    return interp_xs, interp_unc_xs



#_________________________Energy info from another file________________________
foil_energy_data = {'Ni01': {'energy': 27.337003000000003, 'min_unc': 0.6070030000000024, 'plus_unc': 0.6529969999999956}, 
                    'Zr01': {'energy': 26.375581480000005, 'min_unc': 0.605581480000005, 'plus_unc': 0.5944185199999943}, 
                    'Ti01': {'energy': 25.518248420000003, 'min_unc': 0.5882484200000029, 'plus_unc': 0.6117515799999964}, 
                    'Ni02': {'energy': 20.392268865, 'min_unc': 0.7422688649999998, 'plus_unc': 0.7577311350000002}, 
                    'Zr02': {'energy': 19.1745053, 'min_unc': 0.7845052999999993, 'plus_unc': 0.7754946999999994}, 
                    'Ti02': {'energy': 18.058134579999997, 'min_unc': 0.8081345799999973, 'plus_unc': 0.8118654200000037}, 
                    'Ni03': {'energy': 13.9680305, 'min_unc': 0.9780305000000009, 'plus_unc': 1.0619695}, 
                    'Zr03': {'energy': 12.309193775999999, 'min_unc': 1.059193775999999, 'plus_unc': 1.1008062240000012}, 
                    'Ti03': {'energy': 10.706294196000002, 'min_unc': 1.1362941960000015, 'plus_unc': 1.2037058039999984}, 
                    'Ni04': {'energy': 8.648899379944393, 'min_unc': 1.2988993799443929, 'plus_unc': 1.5211006200556056}, 
                    'Zr04': {'energy': 6.082350534817555, 'min_unc': 1.4923505348175548, 'plus_unc': 1.8676494651824447}, 
                    'Ti04': {'energy': 3.6998866038357203, 'min_unc': 2.2298866038357206, 'plus_unc': 2.090113396164279}, 
                    'Ni05': {'energy': 2.2266116583208797, 'min_unc': 1.8966116583208799, 'plus_unc': 1.4033883416791202}, 
                    'Zr05': {'energy': 1.5675066930508672, 'min_unc': 1.2375066930508674, 'plus_unc': 0.7424933069491328}, 
                    'Ti05': {'energy': 0.7071428571428571, 'min_unc': 0.37714285714285717, 'plus_unc': 0.2828571428571429}}




#__________________________Running the function________________________________

# dp = 0.950
# calc_xs_list_of_list_Ni, calc_xs_unc_list_of_list_Ni, beam_energy_in_foil_list_list_Ni, reaction_list_list_Ni = caclulate_xs_in_foil(dp, 'Ni')
# calc_xs_list_of_list_Ti, calc_xs_unc_list_of_list_Ti, beam_energy_in_foil_list_list_Ti, reaction_list_list_Ti = caclulate_xs_in_foil(dp, 'Ti')

# print('_____________p0, dp= 0.988___________________')
calc_xs_list_of_list_Ni_p0, calc_xs_unc_list_of_list_Ni_p0, beam_energy_in_foil_list_list_Ni_p0, reaction_list_list_Ni_p0 = caclulate_xs_in_foil(0.986, 'Ni')
calc_xs_list_of_list_Ti_p0, calc_xs_unc_list_of_list_Ti_p0, beam_energy_in_foil_list_list_Ti_p0, reaction_list_list_Ti_p0 = caclulate_xs_in_foil(0.986, 'Ti')
# print('_____________p1, dp= 0.977___________________')
calc_xs_list_of_list_Ni_p1, calc_xs_unc_list_of_list_Ni_p1, beam_energy_in_foil_list_list_Ni_p1, reaction_list_list_Ni_p1 = caclulate_xs_in_foil(0.972, 'Ni')
calc_xs_list_of_list_Ti_p1, calc_xs_unc_list_of_list_Ti_p1, beam_energy_in_foil_list_list_Ti_p1, reaction_list_list_Ti_p1 = caclulate_xs_in_foil(0.972, 'Ti')




#___________________________Sorting the data___________________________________

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


data_by_reaction_p0 = {}

# Iterate through reaction_list_list_Ni
for i, reaction_list in enumerate(reaction_list_list_Ni_p0):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p0:
            data_by_reaction_p0[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p0[reaction]['calc_xs'].append(calc_xs_list_of_list_Ni_p0[i][j])
        data_by_reaction_p0[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ni_p0[i][j])
        data_by_reaction_p0[reaction]['energy'].append(beam_energy_in_foil_list_list_Ni_p0[i][j])



# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti_p0):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p0:
            data_by_reaction_p0[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p0[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_p0[i][j])
        data_by_reaction_p0[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_p0[i][j])
        data_by_reaction_p0[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_p0[i][j])

#New values after the scaling factor is applied 
data_by_reaction_p0['58CO']['calc_xs'][-1] = 1.0192801735000916
data_by_reaction_p0['58CO']['calc_xs_unc'][-1] = 0.09889413863949252
data_by_reaction_p0['61CU']['calc_xs'][-1] = 0.03757637008354755
data_by_reaction_p0['61CU']['calc_xs_unc'][-1] = 0.009638847896506993



data_by_reaction_p1 = {}

# Iterate through reaction_list_list_Ni
for i, reaction_list in enumerate(reaction_list_list_Ni_p1):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p1:
            data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Ni_p1[i][j])
        data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ni_p1[i][j])
        data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Ni_p1[i][j])



# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti_p1):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p1:
            data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_p1[i][j])
        data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_p1[i][j])
        data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_p1[i][j])


#New values after the scaling factor is applied 
data_by_reaction_p1['58CO']['calc_xs'][-1] = 1.3668111126405154
data_by_reaction_p1['58CO']['calc_xs_unc'][-1] = 0.13261280969184164
data_by_reaction_p1['61CU']['calc_xs'][-1] = 0.05038830494124271
data_by_reaction_p1['61CU']['calc_xs_unc'][-1] = 0.012925282724536654


#____________________________Getting monitor xs__________________________________
E_mon_array = np.linspace(0,30,100000)

# , xs_list_46SC, xs_unc_list_46SC = get_IAEA_monitro_xs('46SC')
#  xs_list_48V, xs_unc_list_48V = get_IAEA_monitro_xs('48V')
# , xs_list_56CO, xs_unc_list_56CO = get_IAEA_monitro_xs('56CO')
# , xs_list_58CO, xs_unc_list_58CO = get_IAEA_monitro_xs('58CO')
# , xs_list_61CU, xs_unc_list_61CU = get_IAEA_monitro_xs('61CU')
interp_xs_46SC, interp_unc_xs_46SC = get_IAEA_monitro_xs('46SC')
interp_xs_48V, interp_unc_xs_48V = get_IAEA_monitro_xs('48V')
interp_xs_56CO, interp_unc_xs_56CO = get_IAEA_monitro_xs('56CO')
interp_xs_58CO, interp_unc_xs_58CO = get_IAEA_monitro_xs('58CO')
interp_xs_61CU, interp_unc_xs_61CU = get_IAEA_monitro_xs('61CU')






# ______________________________Plotting________________________________________

# marker_list = ['d', '*', 's', '<', 'o']
# color_list = ['deepskyblue', 'mediumseagreen', 'gold', 'violet', 'mediumvioletred']


# # Loop over every reaction
# for i, (reaction, data) in enumerate(data_by_reaction.items()):
#     calc_xs = data['calc_xs']
#     calc_xs_unc = data['calc_xs_unc']
#     energies = data['energy']
#     plt.errorbar(energies, calc_xs, yerr=calc_xs_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=f'{reaction} calculated with avrg bc')

# plt.plot(, xs_list_56CO, color='deepskyblue', label='IAEA 56Co')
# plt.plot(, xs_list_58CO, color='mediumseagreen', label='IAEA 58Co')
# plt.plot(, xs_list_61CU, color='gold', label='IAEA 61Cu')
# plt.plot( xs_list_48V, color='violet', label='IAEA 48V')
# plt.plot(, xs_list_46SC, color='mediumvioletred', label='IAEA 46Sc')

# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'dp = {dp:.3f}')
# plt.legend()
# plt.show()

Ni_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Ni')]
Ni_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Ni')]
Ti_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Ti')]
Ti_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Ti')]


calc_xs_56CO_p0 = data_by_reaction_p0['56CO']['calc_xs']
calc_xs_unc_56CO_p0 = data_by_reaction_p0['56CO']['calc_xs_unc']
energies_56CO_p0 = data_by_reaction_p0['56CO']['energy']

calc_xs_58CO_p0 = data_by_reaction_p0['58CO']['calc_xs']
calc_xs_unc_58CO_p0 = data_by_reaction_p0['58CO']['calc_xs_unc']
energies_58CO_p0 = data_by_reaction_p0['58CO']['energy']

calc_xs_61CU_p0 = data_by_reaction_p0['61CU']['calc_xs']
calc_xs_unc_61CU_p0 = data_by_reaction_p0['61CU']['calc_xs_unc']
energies_61CU_p0 = data_by_reaction_p0['61CU']['energy']

calc_xs_46SC_p0 = data_by_reaction_p0['46SC']['calc_xs']
calc_xs_unc_46SC_p0 = data_by_reaction_p0['46SC']['calc_xs_unc']
energies_46SC_p0 = data_by_reaction_p0['46SC']['energy']

calc_xs_48V_p0 = data_by_reaction_p0['48V']['calc_xs']
calc_xs_unc_48V_p0 = data_by_reaction_p0['48V']['calc_xs_unc']
energies_48V_p0 = data_by_reaction_p0['48V']['energy']


calc_xs_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs']
calc_xs_unc_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs_unc']
energies_56CO_p1 = data_by_reaction_p1['56CO']['energy']

calc_xs_58CO_p1 = data_by_reaction_p1['58CO']['calc_xs']
calc_xs_unc_58CO_p1 = data_by_reaction_p1['58CO']['calc_xs_unc']
energies_58CO_p1 = data_by_reaction_p1['58CO']['energy']

calc_xs_61CU_p1 = data_by_reaction_p1['61CU']['calc_xs']
calc_xs_unc_61CU_p1 = data_by_reaction_p1['61CU']['calc_xs_unc']
energies_61CU_p1 = data_by_reaction_p1['61CU']['energy']

calc_xs_46SC_p1 = data_by_reaction_p1['46SC']['calc_xs']
calc_xs_unc_46SC_p1 = data_by_reaction_p1['46SC']['calc_xs_unc']
energies_46SC_p1 = data_by_reaction_p1['46SC']['energy']

calc_xs_48V_p1 = data_by_reaction_p1['48V']['calc_xs']
calc_xs_unc_48V_p1 = data_by_reaction_p1['48V']['calc_xs_unc']
energies_48V_p1 = data_by_reaction_p1['48V']['energy']


plt.figure(figsize=(8, 6))
plt.errorbar(energies_56CO_p0, calc_xs_56CO_p0, xerr=[Ni_energy_min_unc_list[:-1], Ni_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_56CO_p1, calc_xs_56CO_p1, xerr=[Ni_energy_min_unc_list[:-1], Ni_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_array, interp_xs_56CO(E_mon_array), color='hotpink', label='IAEA')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{56}$Co', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.savefig('./Figures/xs_mon_30MeV_56Co.pdf', dpi=600)
plt.show()

plt.figure(figsize=(8, 6))
plt.errorbar(energies_58CO_p0, calc_xs_58CO_p0, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_58CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_58CO_p1, calc_xs_58CO_p1, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_58CO_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_array, interp_xs_58CO(E_mon_array), color='hotpink', label='IAEA')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{58}$Co', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.savefig('./Figures/xs_mon_30MeV_58Co.pdf', dpi=600)
plt.show()

plt.figure(figsize=(8, 6))
plt.errorbar(energies_61CU_p0, calc_xs_61CU_p0, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_61CU_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_61CU_p1, calc_xs_61CU_p1, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_61CU_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_array, interp_xs_61CU(E_mon_array), color='hotpink', label='IAEA')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{61}$Cu', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.savefig('./Figures/xs_mon_30MeV_61Cu.pdf', dpi=600)
plt.show()

plt.figure(figsize=(8, 6))
plt.errorbar(energies_46SC_p0, calc_xs_46SC_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_46SC_p1, calc_xs_46SC_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_array, interp_xs_46SC(E_mon_array), color='hotpink', label='IAEA')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ti(d,x)$^{46}$Sc', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.savefig('./Figures/xs_mon_30MeV_46Sc.pdf', dpi=600)
plt.show()

plt.figure(figsize=(8, 6))
plt.errorbar(energies_48V_p0, calc_xs_48V_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_48V_p1, calc_xs_48V_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_array, interp_xs_48V(E_mon_array), color='hotpink', label='IAEA')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ti(d,x)$^{48}$V', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.savefig('./Figures/xs_mon_30MeV_48V.pdf', dpi=600)
plt.show()




