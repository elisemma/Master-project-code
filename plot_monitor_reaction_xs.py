import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 
from scipy.interpolate import PchipInterpolator
import re


"""
Okey så jeg skjønner ikke hva som skjer. Etter at jeg la inn data by reaction for hånd stemmer ikke dimensjoner og sånn. Må prøve å gå tilbake i morgen og finne ut av hvor ting går galt. Kanskje lage et nytt plotte scrtipt som ikke henter ting i foil class men får alt fra denne og 50 MeV fila (en vendlig manuell måte). Første jeg kan gjøre er å kommentere ut data by reaction som er limt inn for hånd for å se om det hjelper. og bare jobbe meg tilb ake til når ting fungerte
"""




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
        print('target_material', target_material)
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






def get_IAEA_monitro_xs(target_material, reaction_product):
    filename = './Monitor_cross_section_data/IAEA_monitor_xs_' + target_material + '_dx_' + reaction_product + '.txt'
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
#                     'Ti05': {'energy': 0.7071428571428571, 'min_unc': 0.37714285714285717, 'plus_unc': 0.2828571428571429},
#                     'Fe01': {'energy': 48.25110370000002, 'min_unc': 0.9011037000000144, 'plus_unc': 0.8988962999999828}, 
#                     'Zr06': {'energy': 47.66121943333334, 'min_unc': 0.9112194333333434, 'plus_unc': 0.8887805666666537}, 
#                     'Ti06': {'energy': 47.13022093333334, 'min_unc': 0.8802209333333408, 'plus_unc': 0.9197790666666563}, 
#                     'Fe02': {'energy': 41.876494533333336, 'min_unc': 1.0264945333333273, 'plus_unc': 1.073505466666667}, 
#                     'Zr07': {'energy': 41.21671893333334, 'min_unc': 1.066718933333334, 'plus_unc': 1.0332810666666603}, 
#                     'Ti08': {'energy': 40.62009246666667, 'min_unc': 1.0700924666666722, 'plus_unc': 1.0299075333333363}, 
#                     'Fe03': {'energy': 37.24984646666667, 'min_unc': 1.099846466666662, 'plus_unc': 1.1001535333333408}, 
#                     'Zr08': {'energy': 36.523388133333334, 'min_unc': 1.173388133333333, 'plus_unc': 1.1266118666666713}, 
#                     'Ti09': {'energy': 35.86358033333333, 'min_unc': 1.1135803333333314, 'plus_unc': 1.1864196666666658}, 
#                     'Fe04': {'energy': 32.115906233333334, 'min_unc': 1.2659062333333324, 'plus_unc': 1.2340937666666676}, 
#                     'Zr09': {'energy': 31.29907423333333, 'min_unc': 1.249074233333328, 'plus_unc': 1.2509257666666684}, 
#                     'Ti10': {'energy': 30.5522245, 'min_unc': 1.3022245000000012, 'plus_unc': 1.2977755000000002}, 
#                     'Fe05': {'energy': 26.2370075, 'min_unc': 1.4870075000000007, 'plus_unc': 1.5129924999999993}, 
#                     'Zr10': {'energy': 25.275796966666668, 'min_unc': 1.5257969666666682, 'plus_unc': 1.5742030333333332}, 
#                     'Ti11': {'energy': 24.390982466666664, 'min_unc': 1.5409824666666623, 'plus_unc': 1.5590175333333391}}

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
                    # 'Fe01': {'energy': 48.25110370000002, 'min_unc': 0.9011037000000144, 'plus_unc': 0.8988962999999828}, 
                    # 'Zr06': {'energy': 47.66121943333334, 'min_unc': 0.9112194333333434, 'plus_unc': 0.8887805666666537}, 
                    # 'Ti06': {'energy': 47.13022093333334, 'min_unc': 0.8802209333333408, 'plus_unc': 0.9197790666666563}, 
                    # 'Fe02': {'energy': 41.876494533333336, 'min_unc': 1.0264945333333273, 'plus_unc': 1.073505466666667}, 
                    # 'Zr07': {'energy': 41.21671893333334, 'min_unc': 1.066718933333334, 'plus_unc': 1.0332810666666603}, 
                    # 'Ti08': {'energy': 40.62009246666667, 'min_unc': 1.0700924666666722, 'plus_unc': 1.0299075333333363}, 
                    # 'Fe03': {'energy': 37.24984646666667, 'min_unc': 1.099846466666662, 'plus_unc': 1.1001535333333408}, 
                    # 'Zr08': {'energy': 36.523388133333334, 'min_unc': 1.173388133333333, 'plus_unc': 1.1266118666666713}, 
                    # 'Ti09': {'energy': 35.86358033333333, 'min_unc': 1.1135803333333314, 'plus_unc': 1.1864196666666658}, 
                    # 'Fe04': {'energy': 32.115906233333334, 'min_unc': 1.2659062333333324, 'plus_unc': 1.2340937666666676}, 
                    # 'Zr09': {'energy': 31.29907423333333, 'min_unc': 1.249074233333328, 'plus_unc': 1.2509257666666684}, 
                    # 'Ti10': {'energy': 30.5522245, 'min_unc': 1.3022245000000012, 'plus_unc': 1.2977755000000002}, 
                    # 'Fe05': {'energy': 26.2370075, 'min_unc': 1.4870075000000007, 'plus_unc': 1.5129924999999993}, 
                    # 'Zr10': {'energy': 25.275796966666668, 'min_unc': 1.5257969666666682, 'plus_unc': 1.5742030333333332}, 
                    # 'Ti11': {'energy': 24.390982466666664, 'min_unc': 1.5409824666666623, 'plus_unc': 1.5590175333333391}}




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
# calc_xs_list_of_list_Fe_p1, calc_xs_unc_list_of_list_Fe_p1, beam_energy_in_foil_list_list_Fe_p1, reaction_list_list_Fe_p1 = caclulate_xs_in_foil(0.990, 'Fe')
# calc_xs_list_of_list_Ti_p1_50, calc_xs_unc_list_of_list_Ti_p1_50, beam_energy_in_foil_list_list_Ti_p1_50, reaction_list_list_Ti_p1_50 = caclulate_xs_in_foil(0.990, 'Ti')




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

# # Iterate through reaction_list_list_Fe
# for i, reaction_list in enumerate(reaction_list_list_Fe_p1):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_p1:
#             data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Fe_p1[i][j])
#         data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Fe_p1[i][j])
#         data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Fe_p1[i][j])



# # Iterate through reaction_list_list_Ti 50 MeV
# for i, reaction_list in enumerate(reaction_list_list_Ti_p1_50):
#     for j, reaction in enumerate(reaction_list):
#         if reaction not in data_by_reaction_p1:
#             data_by_reaction_p1[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
#         # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
#         data_by_reaction_p1[reaction]['calc_xs'].append(calc_xs_list_of_list_Ti_p1_50[i][j])
#         data_by_reaction_p1[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list_Ti_p1_50[i][j])
#         data_by_reaction_p1[reaction]['energy'].append(beam_energy_in_foil_list_list_Ti_p1_50[i][j])



#New values after the scaling factor is applied 
data_by_reaction_p1['58CO']['calc_xs'][-1] = 1.3668111126405154
data_by_reaction_p1['58CO']['calc_xs_unc'][-1] = 0.13261280969184164
data_by_reaction_p1['61CU']['calc_xs'][-1] = 0.05038830494124271
data_by_reaction_p1['61CU']['calc_xs_unc'][-1] = 0.012925282724536654

print('_____________________________________')
print(data_by_reaction_p1)
print('_____________________________________')


# data_by_reaction_p1 = {'56CO': {'calc_xs': [np.float64(10.534973121443173), np.float64(8.468968382445144), np.float64(21.40976480170361), np.float64(30.456223871264914)],             
#                                 'calc_xs_unc': [np.float64(2.023700489023643), np.float64(1.4182512775939868), np.float64(1.8112243925372662), np.float64(8.256030363968456)],              
#                                 'energy': [np.float64(27.31117938), np.float64(20.287392734999997), np.float64(13.757549620000002), np.float64(8.280190492265563)]}, 
#                        '58CO': {'calc_xs': [np.float64(137.91150626639745), np.float64(226.96776626771035), np.float64(180.63987499745105), np.float64(47.589777075292574), 1.0192801735000916], 
#                                 'calc_xs_unc': [np.float64(25.964617540369627), np.float64(50.350128742727186), np.float64(14.96681581994496), np.float64(21.703710247522483), 0.09889413863949252], 
#                                 'energy': [np.float64(27.31117938), np.float64(20.287392734999997), np.float64(13.757549620000002), np.float64(8.280190492265563), np.float64(2.0589693651507943)]}, 
#                        '61CU': {'calc_xs': [np.float64(15.424367205216003), np.float64(15.887056603664542), np.float64(25.695573663908974), np.float64(57.76418305132154), 0.03757637008354755], 
#                                 'calc_xs_unc': [np.float64(2.7941276596308304), np.float64(2.6753550655162144), np.float64(2.171368266450643), np.float64(15.689424358708196), 0.009638847896506993], 
#                                 'energy': [np.float64(27.31117938), np.float64(20.287392734999997), np.float64(13.757549620000002), np.float64(8.280190492265563), np.float64(2.0589693651507943)]}, 
                       
#                        '48V': {'calc_xs': [np.float64(53.28131917888242), np.float64(71.53510820381517), np.float64(86.95500752412562), np.float64(120.32118113088468), np.float64(211.09294415195407), np.float64(186.9372369978318), np.float64(328.81764193119716), np.float64(140.34536460253315), np.float64(10.611539466199133)], 
#                                'calc_xs_unc': [np.float64(4.004285252699319), np.float64(2.215160347196869), np.float64(2.3840605321553276), np.float64(5.720374572911209), np.float64(2.7645225694942264), np.float64(6.451035557195398), np.float64(13.197687966450168), np.float64(11.493864213759096), np.float64(4.927076297067808)], 
#                                'energy': [np.float64(47.09795026666667), np.float64(40.50592943333334), np.float64(35.68120613333333), np.float64(30.28116373333334), np.float64(23.99016156666667), np.float64(25.473729640000002), np.float64(17.92027038), np.float64(10.413339984), np.float64(3.3300379128898365)]}, 
                       
#                        '46SC': {'calc_xs': [np.float64(71.23392223363822), np.float64(47.64071433003223), np.float64(31.38176508759792), np.float64(25.349146020414114), np.float64(23.97580225065507), np.float64(23.90542380806754), np.float64(25.50009762253962), np.float64(30.90276823909987), np.float64(1.9078634847081497)], 
#                                 'calc_xs_unc': [np.float64(5.412575994612804), np.float64(1.4596181121256113), np.float64(0.8335479522980904), np.float64(1.1962400980743435), np.float64(0.3986503440442611), np.float64(0.8545507969898345), np.float64(1.0256488457715143), np.float64(2.5372147352448016), np.float64(0.8874467747252316)], 
#                                 'energy': [np.float64(47.09795026666667), np.float64(40.50592943333334), np.float64(35.68120613333333), np.float64(30.28116373333334), np.float64(23.99016156666667), np.float64(25.473729640000002), np.float64(17.92027038), np.float64(10.413339984), np.float64(3.3300379128898365)]}}

# {'56CO': {'calc_xs': [np.float64(42.80199907072412), np.float64(53.53658424185413), np.float64(67.42712662562114), np.float64(95.89102219731262), np.float64(165.9191705527795)], 'calc_xs_unc': [np.float64(3.0335331140655164), np.float64(3.843651004670334), np.float64(4.898118323109198), np.float64(6.741815080623053), np.float64(11.665288091765959)], 'energy': [np.float64(48.2318781), np.float64(41.779047799999994), np.float64(37.08828046666665), np.float64(31.872687066666668), np.float64(25.879117466666663)]}, 
 
#  '46SC': {'calc_xs': [np.float64(71.23392223363822), np.float64(47.64071433003223), np.float64(31.38176508759792), np.float64(25.349146020414114), np.float64(23.97580225065507)], 
#           'calc_xs_unc': [np.float64(5.412575994612804), np.float64(1.4596181121256113), np.float64(0.8335479522980904), np.float64(1.1962400980743435), np.float64(0.3986503440442611)],
#             'energy': [np.float64(47.09795026666667), np.float64(40.50592943333334), np.float64(35.68120613333333), np.float64(30.28116373333334), np.float64(23.99016156666667)]}, 
 
# '48V': {'calc_xs': [np.float64(53.28131917888242), np.float64(71.53510820381517), np.float64(86.95500752412562), np.float64(120.32118113088468), np.float64(211.09294415195407)], 
#          'calc_xs_unc': [np.float64(4.004285252699319), np.float64(2.215160347196869), np.float64(2.3840605321553276), np.float64(5.720374572911209), np.float64(2.7645225694942264)], 
#          'energy': [np.float64(47.09795026666667), np.float64(40.50592943333334), np.float64(35.68120613333333), np.float64(30.28116373333334), np.float64(23.99016156666667)]}}


#____________________________Getting monitor xs__________________________________
E_mon_array = np.linspace(0,50,100000)

# , xs_list_46SC, xs_unc_list_46SC = get_IAEA_monitro_xs('46SC')
#  xs_list_48V, xs_unc_list_48V = get_IAEA_monitro_xs('48V')
# , xs_list_56CO, xs_unc_list_56CO = get_IAEA_monitro_xs('56CO')
# , xs_list_58CO, xs_unc_list_58CO = get_IAEA_monitro_xs('58CO')
# , xs_list_61CU, xs_unc_list_61CU = get_IAEA_monitro_xs('61CU')
interp_xs_Ti_dx_46SC, interp_unc_xs_Ti_dx_46SC = get_IAEA_monitro_xs('Ti','46SC')
interp_xs_Ti_dx_48V, interp_unc_xs_Ti_dx_48V = get_IAEA_monitro_xs('Ti','48V')
interp_xs_Ni_dx_56CO, interp_unc_xs_Ni_dx_56CO = get_IAEA_monitro_xs('Ni','56CO')
interp_xs_Ni_dx_58CO, interp_unc_xs_Ni_dx_58CO = get_IAEA_monitro_xs('Ni','58CO')
interp_xs_Ni_dx_61CU, interp_unc_xs_Ni_dx_61CU = get_IAEA_monitro_xs('Ni','61CU')
# interp_xs_Fe_dx_56CO, interp_unc_xs_Fe_dx_56CO = get_IAEA_monitro_xs('Fe','56CO')






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
Fe_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Fe')]
Fe_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Fe')]

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

# calc_xs_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs']
# calc_xs_unc_56CO_p1 = data_by_reaction_p1['56CO']['calc_xs_unc'] #Må skille på Fe og Ni!! xxx
# energies_56CO_p1 = data_by_reaction_p1['56CO']['energy']



#Get data from exfor
markers = ['.', '*', 'v', '^', '+', '<', '>', 's', 'h',     '.', '*', 'v', '^', '+', '<', '>', 's', 'h']
grey_colors = ['dimgrey', 'darkgrey', 'lightgrey', 'silver', 'k', 'dimgrey', 'darkgrey', 'lightgrey', 'silver',     'silver', 'lightgrey', 'darkgrey', 'dimgrey', 'k','silver', 'lightgrey', 'dimgrey', 'darkgrey']









plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ni_dx_56CO/'
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year})"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_56CO_p0, calc_xs_56CO_p0, xerr=[Ni_energy_min_unc_list[:-1], Ni_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ni_dx_56CO(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ni_dx_56CO(E_mon_array)-interp_unc_xs_Ni_dx_56CO(E_mon_array), interp_xs_Ni_dx_56CO(E_mon_array)+interp_unc_xs_Ni_dx_56CO(E_mon_array), color='lavender')
plt.errorbar(energies_56CO_p1, calc_xs_56CO_p1, xerr=[Ni_energy_min_unc_list[:-1], Ni_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{56}$Co', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.ylim(0,40)
# plt.savefig('./Figures/xs_mon_30MeV_56Co.pdf', dpi=600)
plt.show()







plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ni_dx_58CO/'
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year})"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_58CO_p0, calc_xs_58CO_p0, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_58CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ni_dx_58CO(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ni_dx_58CO(E_mon_array)-interp_unc_xs_Ni_dx_58CO(E_mon_array), interp_xs_Ni_dx_58CO(E_mon_array)+interp_unc_xs_Ni_dx_58CO(E_mon_array), color='lavender')
plt.errorbar(energies_58CO_p1, calc_xs_58CO_p1, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_58CO_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{58}$Co', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.ylim(0)
# plt.savefig('./Figures/xs_mon_30MeV_58Co.pdf', dpi=600)
plt.show()







plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ni_dx_61CU/'
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year})"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_61CU_p0, calc_xs_61CU_p0, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_61CU_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ni_dx_61CU(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ni_dx_61CU(E_mon_array)-interp_unc_xs_Ni_dx_61CU(E_mon_array), interp_xs_Ni_dx_61CU(E_mon_array)+interp_unc_xs_Ni_dx_61CU(E_mon_array), color='lavender')
plt.errorbar(energies_61CU_p1, calc_xs_61CU_p1, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_61CU_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{61}$Cu', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.ylim(0)
# plt.savefig('./Figures/xs_mon_30MeV_61Cu.pdf', dpi=600)
plt.show()







plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ti_dx_46SC/'
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year})"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_46SC_p0, calc_xs_46SC_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ti_dx_46SC(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ti_dx_46SC(E_mon_array)-interp_unc_xs_Ti_dx_46SC(E_mon_array), interp_xs_Ti_dx_46SC(E_mon_array)+interp_unc_xs_Ti_dx_46SC(E_mon_array), color='lavender')
plt.errorbar(energies_46SC_p1, calc_xs_46SC_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ti(d,x)$^{46}$Sc', fontsize=14)
plt.legend()
plt.xlim(0,50)
plt.ylim(0)
# plt.savefig('./Figures/xs_mon_30MeV_46Sc.pdf', dpi=600)
plt.show()






plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ti_dx_48V/'
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year})"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_48V_p0, calc_xs_48V_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ti_dx_48V(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ti_dx_48V(E_mon_array)-interp_unc_xs_Ti_dx_48V(E_mon_array), interp_xs_Ti_dx_48V(E_mon_array)+interp_unc_xs_Ti_dx_48V(E_mon_array), color='lavender')
plt.errorbar(energies_48V_p1, calc_xs_48V_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ti(d,x)$^{48}$V', fontsize=14)
plt.legend()
plt.xlim(0,50)
plt.ylim(0)
# plt.savefig('./Figures/xs_mon_30MeV_48V.pdf', dpi=600)
plt.show()





# plt.figure(figsize=(8, 6))
# path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Fe_dx_56CO/'
#  # List all files in the folder
# files = os.listdir(path)
# # Loop through each file
# for i, file_name in enumerate(files):
#     file_path = os.path.join(path, file_name)
#     data = np.loadtxt(file_path, comments=["#", "//"])
#     energy = data[:, 0]
#     energy_unc = data[:, 1]
#     cross_section = data[:, 2] * 1e3
#     cross_section_unc = data[:, 3] * 1e3
#     # Initialize variables to store author name and year
#     author_name = ""
#     year = ""
#     # Read the file to find author name and year
#     with open(file_path, 'r') as file:
#         for line in file:
#             if '# ' in line:  # Check if the line contains a comment
#                 # Use regular expression to find the author name and year
#                 match = re.search(r'(\d{4}),([^#]+)', line)
#                 if match:
#                     year = match.group(1)  # Extract the year
#                     author_name = match.group(2).strip()  # Extract the author name
#                     # Remove unwanted single-letter followed by a dot and the plus sign
#                     author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
#                     break  # Stop after finding the first occurrence of name and year
#     # Create the label
#     label = f"{author_name} ({year})"
#     plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.plot(E_mon_array, interp_xs_Fe_dx_56CO(E_mon_array), color='cornflowerblue', label='IAEA')
# plt.fill_between(E_mon_array, interp_xs_Fe_dx_56CO(E_mon_array)-interp_unc_xs_Fe_dx_56CO(E_mon_array), interp_xs_Fe_dx_56CO(E_mon_array)+interp_unc_xs_Fe_dx_56CO(E_mon_array), color='lavender')
# plt.errorbar(energies_56CO_p1, calc_xs_56CO_p1, xerr=[Fe_energy_min_unc_list[:-1], Fe_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
# plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
# plt.ylabel('Cross Section (mb)', fontsize=14)
# plt.title(r'$^{nat}$Fe(d,x)$^{56}$Co', fontsize=14)
# plt.legend()
# plt.xlim(0,50)
# plt.ylim(0)
# # plt.savefig('./Figures/xs_mon_30MeV_56Co.pdf', dpi=600)
# plt.show()



