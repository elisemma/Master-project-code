import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 



def caclulate_xs_in_foil(dp, compound):
    stack_df = pd.read_csv(f'./Stack_calculations/stack_50MeV_dp_{dp:.3f}.csv')
    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]

    calc_xs_list_of_list = [] 
    calc_xs_unc_list_of_list = []

    beam_energy_in_foil_list_list = []

    reaction_list_list = []


    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']
        # beam_energy_in_foil = row['mu_E']

        A0_by_curie_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_curie.csv')
        # if os.path.exists(f'./Calculated_A0/{foil_name}_A0_by_hand.csv'):
        #     A0_by_hand_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_hand.csv')
        #     A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)
        # else:
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



#_________________________Energy info from another file________________________
foil_energy_data ={'Fe01': {'energy': 48.25110370000002, 'min_unc': 0.9011037000000144, 'plus_unc': 0.8988962999999828}, 
                   'Zr06': {'energy': 47.66121943333334, 'min_unc': 0.9112194333333434, 'plus_unc': 0.8887805666666537}, 
                   'Ti06': {'energy': 47.13022093333334, 'min_unc': 0.8802209333333408, 'plus_unc': 0.9197790666666563}, 
                   'Fe02': {'energy': 41.876494533333336, 'min_unc': 1.0264945333333273, 'plus_unc': 1.073505466666667}, 
                   'Zr07': {'energy': 41.21671893333334, 'min_unc': 1.066718933333334, 'plus_unc': 1.0332810666666603}, 
                   'Ti08': {'energy': 40.62009246666667, 'min_unc': 1.0700924666666722, 'plus_unc': 1.0299075333333363}, 
                   'Fe03': {'energy': 37.24984646666667, 'min_unc': 1.099846466666662, 'plus_unc': 1.1001535333333408}, 
                   'Zr08': {'energy': 36.523388133333334, 'min_unc': 1.173388133333333, 'plus_unc': 1.1266118666666713}, 
                   'Ti09': {'energy': 35.86358033333333, 'min_unc': 1.1135803333333314, 'plus_unc': 1.1864196666666658}, 
                   'Fe04': {'energy': 32.115906233333334, 'min_unc': 1.2659062333333324, 'plus_unc': 1.2340937666666676}, 
                   'Zr09': {'energy': 31.29907423333333, 'min_unc': 1.249074233333328, 'plus_unc': 1.2509257666666684}, 
                   'Ti10': {'energy': 30.5522245, 'min_unc': 1.3022245000000012, 'plus_unc': 1.2977755000000002}, 
                   'Fe05': {'energy': 26.2370075, 'min_unc': 1.4870075000000007, 'plus_unc': 1.5129924999999993}, 
                   'Zr10': {'energy': 25.275796966666668, 'min_unc': 1.5257969666666682, 'plus_unc': 1.5742030333333332}, 
                   'Ti11': {'energy': 24.390982466666664, 'min_unc': 1.5409824666666623, 'plus_unc': 1.5590175333333391}}




#__________________________Running the function________________________________

# dp = 0.97
# calc_xs_list_of_list_Fe, calc_xs_unc_list_of_list_Fe, beam_energy_in_foil_list_list_Fe, reaction_list_list_Fe = caclulate_xs_in_foil(dp, 'Fe')
# calc_xs_list_of_list_Ti, calc_xs_unc_list_of_list_Ti, beam_energy_in_foil_list_list_Ti, reaction_list_list_Ti = caclulate_xs_in_foil(dp, 'Ti')


calc_xs_list_of_list_Fe_p0, calc_xs_unc_list_of_list_Fe_p0, beam_energy_in_foil_list_list_Fe_p0, reaction_list_list_Fe_p0 = caclulate_xs_in_foil(0.965, 'Fe')
calc_xs_list_of_list_Ti_p0, calc_xs_unc_list_of_list_Ti_p0, beam_energy_in_foil_list_list_Ti_p0, reaction_list_list_Ti_p0 = caclulate_xs_in_foil(0.965, 'Ti')

calc_xs_list_of_list_Fe_p1, calc_xs_unc_list_of_list_Fe_p1, beam_energy_in_foil_list_list_Fe_p1, reaction_list_list_Fe_p1 = caclulate_xs_in_foil(0.979, 'Fe')
calc_xs_list_of_list_Ti_p1, calc_xs_unc_list_of_list_Ti_p1, beam_energy_in_foil_list_list_Ti_p1, reaction_list_list_Ti_p1 = caclulate_xs_in_foil(0.979, 'Ti')

# calc_xs_list_of_list_Fe_avrg, calc_xs_unc_list_of_list_Fe_avrg, beam_energy_in_foil_list_list_Fe_avrg, reaction_list_list_Fe_avrg = caclulate_xs_in_foil(0.972, 'Fe')
# calc_xs_list_of_list_Ti_avrg, calc_xs_unc_list_of_list_Ti_avrg, beam_energy_in_foil_list_list_Ti_avrg, reaction_list_list_Ti_avrg = caclulate_xs_in_foil(0.972, 'Ti')


#___________________________Sorting the data___________________________________

# # Initialize nested dictionary to store beam currents, uncertainties and energies for each reaction
# data_by_reaction = {}

# # Iterate through reaction_list_list_Fe
# for i, reaction_list in enumerate(reaction_list_list_Fe):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction:
#             data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list_Fe[i][j])
#         data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Fe[i][j])
#         data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list_Fe[i][j])



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
for i, reaction_list in enumerate(reaction_list_list_Fe_p0):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p0:
            data_by_reaction_p0[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p0[reaction]['calc_xs'].append(calc_xs_list_of_list_Fe_p0[i][j])
        data_by_reaction_p0[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Fe_p0[i][j])
        data_by_reaction_p0[reaction]['energy'].append(beam_energy_in_foil_list_list_Fe_p0[i][j])



# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti_p0):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p0:
            data_by_reaction_p0[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p0[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_p0[i][j])
        data_by_reaction_p0[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_p0[i][j])
        data_by_reaction_p0[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_p0[i][j])


data_by_reaction_p1 = {}

# Iterate through reaction_list_list_Ni
for i, reaction_list in enumerate(reaction_list_list_Fe_p1):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p1:
            data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Fe_p1[i][j])
        data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Fe_p1[i][j])
        data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Fe_p1[i][j])



# Iterate through reaction_list_list_Ti
for i, reaction_list in enumerate(reaction_list_list_Ti_p1):
    for j, reaction in enumerate(reaction_list):
        if reaction not in data_by_reaction_p1:
            data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
        # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
        data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_p1[i][j])
        data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_p1[i][j])
        data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_p1[i][j])




# data_by_reaction_avrg = {}

# # Iterate through reaction_list_list_Ni
# for i, reaction_list in enumerate(reaction_list_list_Fe_avrg):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_avrg:
#             data_by_reaction_avrg[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_avrg[reaction]['calc_xs'].append(calc_xs_list_of_list_Fe_avrg[i][j])
#         data_by_reaction_avrg[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Fe_avrg[i][j])
#         data_by_reaction_avrg[reaction]['energy'].append(beam_energy_in_foil_list_list_Fe_avrg[i][j])



# # Iterate through reaction_list_list_Ti
# for i, reaction_list in enumerate(reaction_list_list_Ti_avrg):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_avrg:
#             data_by_reaction_avrg[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_avrg[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_avrg[i][j])
#         data_by_reaction_avrg[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_avrg[i][j])
#         data_by_reaction_avrg[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_avrg[i][j])



#____________________________Getting monitor xs__________________________________

E_mon_list_46SC, xs_list_46SC, xs_unc_list_46SC = get_IAEA_monitro_xs('46SC')
E_mon_list_48V, xs_list_48V, xs_unc_list_48V = get_IAEA_monitro_xs('48V')
E_mon_list_56CO, xs_list_56CO, xs_unc_list_56CO = get_IAEA_monitro_xs('56CO')






# # ______________________________Plotting________________________________________

# marker_list = ['d', '*', 's']
# color_list = ['deepskyblue','gold', 'hotpink']


# # Loop over every reaction
# for i, (reaction, data) in enumerate(data_by_reaction.items()):
#     calc_xs = data['calc_xs']
#     calc_xs_unc = data['calc_xs_unc']
#     energies = data['energy']
#     plt.errorbar(energies, calc_xs, yerr=calc_xs_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=f'{reaction} calculated with avrg bc')

# plt.plot(E_mon_list_56CO, xs_list_56CO, color='deepskyblue', label='IAEA 56Co')
# plt.plot(E_mon_list_46SC, xs_list_46SC, color='gold', label='IAEA 46Sc')
# plt.plot(E_mon_list_48V, xs_list_48V, color='hotpink', label='IAEA 48V')


# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'dp = {dp:.3f}')
# plt.legend()
# plt.show()

Fe_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Fe')]
Fe_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Fe')]
Ti_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Ti')]
Ti_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Ti')]


calc_xs_56CO_p0 = data_by_reaction_p0['56CO']['calc_xs']
calc_xs_unc_56CO_p0 = data_by_reaction_p0['56CO']['calc_xs_unc']
energies_56CO_p0 = data_by_reaction_p0['56CO']['energy']

calc_xs_46SC_p0 = data_by_reaction_p0['46SC']['calc_xs']
calc_xs_unc_46SC_p0 = data_by_reaction_p0['46SC']['calc_xs_unc']
energies_46SC_p0 = data_by_reaction_p0['46SC']['energy']

calc_xs_48V_p0 = data_by_reaction_p0['48V']['calc_xs']
calc_xs_unc_48V_p0 = data_by_reaction_p0['48V']['calc_xs_unc']
energies_48V_p0 = data_by_reaction_p0['48V']['energy']


calc_xs_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs']
calc_xs_unc_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs_unc']
energies_56CO_p1 = data_by_reaction_p1['56CO']['energy']

calc_xs_46SC_p1 = data_by_reaction_p1['46SC']['calc_xs']
calc_xs_unc_46SC_p1 = data_by_reaction_p1['46SC']['calc_xs_unc']
energies_46SC_p1 = data_by_reaction_p1['46SC']['energy']

calc_xs_48V_p1 = data_by_reaction_p1['48V']['calc_xs']
calc_xs_unc_48V_p1 = data_by_reaction_p1['48V']['calc_xs_unc']
energies_48V_p1 = data_by_reaction_p1['48V']['energy']

# calc_xs_56CO_avrg = data_by_reaction_avrg['56CO']['calc_xs']
# calc_xs_unc_56CO_avrg = data_by_reaction_avrg['56CO']['calc_xs_unc']
# energies_56CO_avrg = data_by_reaction_avrg['56CO']['energy']

# calc_xs_46SC_avrg = data_by_reaction_avrg['46SC']['calc_xs']
# calc_xs_unc_46SC_avrg = data_by_reaction_avrg['46SC']['calc_xs_unc']
# energies_46SC_avrg = data_by_reaction_avrg['46SC']['energy']

# calc_xs_48V_avrg = data_by_reaction_avrg['48V']['calc_xs']
# calc_xs_unc_48V_avrg = data_by_reaction_avrg['48V']['calc_xs_unc']
# energies_48V_avrg = data_by_reaction_avrg['48V']['energy']



plt.errorbar(energies_56CO_p0, calc_xs_56CO_p0, xerr=[Fe_energy_min_unc_list, Fe_energy_plus_unc_list], yerr=calc_xs_unc_56CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_56CO_p1, calc_xs_56CO_p1, xerr=[Fe_energy_min_unc_list, Fe_energy_plus_unc_list], yerr=calc_xs_unc_56CO_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_list_56CO, xs_list_56CO, color='hotpink', label='IAEA')
plt.xlabel('Beam energy (MeV)')
plt.ylabel('Cross section (mb)')
plt.title(f'56CO')
plt.legend()
plt.xlim(0,50)
plt.show()

plt.errorbar(energies_46SC_p0, calc_xs_46SC_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_46SC_p1, calc_xs_46SC_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_list_46SC, xs_list_46SC, color='hotpink', label='IAEA')
plt.xlabel('Beam energy (MeV)')
plt.ylabel('Cross section (mb)')
plt.title(f'46SC')
plt.legend()
plt.xlim(0,50)
plt.show()

plt.errorbar(energies_48V_p0, calc_xs_48V_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.errorbar(energies_48V_p1, calc_xs_48V_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
plt.plot(E_mon_list_48V, xs_list_48V, color='hotpink', label='IAEA')
plt.xlabel('Beam energy (MeV)')
plt.ylabel('Cross section (mb)')
plt.title(f'48V')
plt.legend()
plt.xlim(0,50)
plt.show()


# plt.errorbar(energies_56CO_p0, calc_xs_56CO_p0, yerr=calc_xs_unc_56CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
# plt.errorbar(energies_56CO_p1, calc_xs_56CO_p1, yerr=calc_xs_unc_56CO_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
# # plt.errorbar(energies_56CO_avrg, calc_xs_56CO_avrg, yerr=calc_xs_unc_56CO_avrg, marker='d', markersize=5, linestyle='', color='grey', label='avrg')
# plt.plot(E_mon_list_56CO, xs_list_56CO, color='hotpink', label='IAEA')
# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'56CO')
# plt.legend()
# plt.xlim(0,50)
# plt.show()

# plt.errorbar(energies_46SC_p0, calc_xs_46SC_p0, yerr=calc_xs_unc_46SC_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
# plt.errorbar(energies_46SC_p1, calc_xs_46SC_p1, yerr=calc_xs_unc_46SC_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
# # plt.errorbar(energies_46SC_avrg, calc_xs_46SC_avrg, yerr=calc_xs_unc_46SC_avrg, marker='d', markersize=5, linestyle='', color='grey', label='avrg')
# plt.plot(E_mon_list_46SC, xs_list_46SC, color='hotpink', label='IAEA')
# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'46SC')
# plt.legend()
# plt.xlim(0,50)
# plt.show()

# plt.errorbar(energies_48V_p0, calc_xs_48V_p0, yerr=calc_xs_unc_48V_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
# plt.errorbar(energies_48V_p1, calc_xs_48V_p1, yerr=calc_xs_unc_48V_p1, marker='d', markersize=5, linestyle='', color='gold', label='p1')
# # plt.errorbar(energies_48V_avrg, calc_xs_48V_avrg, yerr=calc_xs_unc_48V_avrg, marker='d', markersize=5, linestyle='', color='grey', label='avrg')
# plt.plot(E_mon_list_48V, xs_list_48V, color='hotpink', label='IAEA')
# plt.xlabel('Beam energy (MeV)')
# plt.ylabel('Cross section (mb)')
# plt.title(f'48V')
# plt.legend()
# plt.xlim(0,50)
# plt.show()




