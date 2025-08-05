import numpy as np 
import matplotlib.pyplot as plt


def calculate_weighted_average_beam_current(beam_current_list, beam_current_unc_list):
    beam_currents = np.array(beam_current_list)
    beam_currents_unc = np.array(beam_current_unc_list)

    weight_list = []

    for i in range(len(beam_currents)):
        weight = 1 / beam_currents_unc[i]**2
        weight_list.append(weight)

    # Calculate the weighted average beam current
    weighted_average_beam_current = np.average(beam_currents, weights=weight_list)

    # Calculate the variance using the weighted_average_beam_current
    weighted_average_beam_current_array = np.zeros(len(beam_currents))
    weighted_average_beam_current_array.fill(weighted_average_beam_current)
    sigma = np.sqrt(np.average((beam_currents - weighted_average_beam_current_array)**2, weights=weight_list))

    return weighted_average_beam_current, sigma






data_by_reaction = {'56CO': {'beam_current': [142.94357895718744, 140.42980469942898, 137.16906671032277, 115.45740940361617], 
                             'beam_current_unc': [16.1669084292634, 17.49449730400099, 16.161155916519576, 13.859310110622088], 
                             'energy': [27.337003000000003, 20.392268865, 13.9680305, 8.648899379944393]}, 
                    '58CO': {'beam_current': [120.23148693034973, 131.67191999206304, 131.75266364143582, 126.1803133367875, 481.99333497949357], 
                             'beam_current_unc': [13.001483412818153, 14.233075670536211, 15.866540603364134, 20.456362413666625, 111.93272788291921], 
                             'energy': [27.337003000000003, 20.392268865, 13.9680305, 8.648899379944393, 2.2266116583208797]}, 
                    '61CU': {'beam_current': [155.59646977043542, 122.99063987093889, 135.6182994747526, 111.09164503202697, 7.3406068653733865], 
                             'beam_current_unc': [16.658675796237205, 12.754575896483685, 12.965154905748303, 13.097625213556269, 1.6499226209787192], 
                             'energy': [27.337003000000003, 20.392268865, 13.9680305, 8.648899379944393, 2.2266116583208797]}, 
                    '48V': {'beam_current': [124.02956425121741, 120.2778569729406, 127.58646137864346, 141.91860731811067, 993.6251140167863], 
                             'beam_current_unc': [13.022205611149174, 12.580416248317658, 13.888015254093956, 30.341940981048147, 288.69582550692076], 
                             'energy': [25.518248420000003, 18.058134579999997, 10.706294196000002, 3.6998866038357203, 0.7071428571428571]}, 
                    '46SC': {'beam_current': [119.74538210590372, 116.4687543732036, 127.20246828201854, 72.52638945931477, 1439.6284842936825], 
                             'beam_current_unc': [10.095801563466392, 9.781591340171868, 11.017872005073386, 8.933645647750188, 738.0634608403565], 
                             'energy': [25.518248420000003, 18.058134579999997, 10.706294196000002, 3.6998866038357203, 0.7071428571428571]}}


weighted_avrg_bc = [125.94248971622208, 129.10847496934917, 131.37180938909623, 113.73832236122058,]
weighted_avrg_bc_unc = [1.285853751774515, 1.1323612222751478, 0.5845008492443324, 1.56124875976737]

avrg_bc = [128.21543781323695, 123.4575028703829, 131.18779932476846,  96.63429803291486]
avrg_bc_unc = [12.74029303652092, 7.710074388016526, 4.0835044453116245, 23.282870237984515]

mu_energy = [26.3, 19.2, 12.3, 6.1] #These are not right

def calculate_std_per_position(data_by_reaction):
    # Find the maximum length among all beam_current lists
    max_length = max(len(data['beam_current']) for data in data_by_reaction.values())

    # Pad shorter lists with np.nan
    for reaction, data in data_by_reaction.items():
        current_length = len(data['beam_current'])
        if current_length < max_length:
            data['beam_current'] += [np.nan] * (max_length - current_length)

    # Calculate standard deviation for each position in the list
    std_dev_per_position = [np.nanstd([data['beam_current'][i] for data in data_by_reaction.values()]) for i in range(max_length)]

    # Print standard deviation for each position
    for i, std_dev in enumerate(std_dev_per_position):
        print(f"Standard deviation of beam current at position {i+1}: {std_dev}")

# calculate_std_per_position(data_by_reaction)


# data_by_reaction = {'56CO': {'beam_current': [91.11279483618009, 95.48943423201047, 96.67759656591667, 87.90385977930694, 79.76353046504703], 
#                              'beam_current_unc': [17.764379182058054, 17.550901526554565, 16.71724750021088, 13.273019225954833, 9.946307812204912], 
#                              'energy': [48.25110370000002, 41.876494533333336, 37.24984646666667, 32.115906233333334, 26.2370075]}, 
#                     '46SC': {'beam_current': [98.38985492917791, 92.43719344140617, 94.87157275894245, 92.64152715822723, 89.76044811394085], 
#                              'beam_current_unc': [9.06208185316921, 8.211563169758909, 8.188146100794526, 8.695843008929224, 7.418752250751914], 
#                              'energy': [47.13022093333334, 40.62009246666667, 35.86358033333333, 30.5522245, 24.390982466666664]}, 
#                     '48V':  {'beam_current': [88.69779045350535, 89.04935458091553, 92.2968314456004, 90.59958822021765, 89.8492617019322], 
#                              'beam_current_unc': [12.20254718194584, 9.989933425556139, 9.93150127659525, 10.305181149242912, 9.318875229511063], 
#                              'energy': [47.13022093333334, 40.62009246666667, 35.86358033333333, 30.5522245, 24.390982466666664]}}



# weighted_avrg_bc = [90.40077287936163, 90.93641924096387, 94.82509814227382, 91.51235172180517, 89.92515787653733]
# weighted_avrg_bc_unc = [1.1399343628389937, 1.018094891664483, 0.776117718006085, 2.0148678905165527, 0.31903965982430116]

# avrg_bc = [94.39425492213508, 91.58165388637218, 94.18557457328106, 91.01283131007776, 87.24670967193768]
# avrg_bc_unc = [4.498503761693877, 2.1065793327169233, 1.5114371572602803, 1.7980871699537482, 4.36691737673634]

# mu_energy = [47.7, 41.2, 36.5, 30.5, 23.5]



# # Calculate standard deviation for each position in the list
# num_elements = len(next(iter(data_by_reaction.values()))['beam_current'])  # Assuming all lists have the same length
# std_dev_per_position = [np.std([data['beam_current'][i] for data in data_by_reaction.values()]) for i in range(num_elements)]

# # Print standard deviation for each position
# for i, std_dev in enumerate(std_dev_per_position):
#     print(f"Standard deviation of beam current in comp {i+1}: {std_dev}")



# ______________________________Plotting________________________________________

marker_list = ['d', '*', 's', '<', 'o']
color_list = ['deepskyblue', 'mediumseagreen', 'gold', 'violet', 'mediumvioletred']

# Loop over every reaction
for i, (reaction, data) in enumerate(data_by_reaction.items()):
    beam_currents = data['beam_current']
    beam_currents_unc = data['beam_current_unc']
    energies = data['energy']

    plt.errorbar(energies, beam_currents, yerr=beam_currents_unc, marker=marker_list[i], markersize=5, linestyle='', color=color_list[i], label=reaction, capsize = 5)

plt.errorbar(mu_energy, weighted_avrg_bc, yerr = weighted_avrg_bc_unc, marker = 's', markersize=7, linestyle='', capsize = 7, color='black', label='weighted avrg')
plt.errorbar(mu_energy, avrg_bc, yerr = avrg_bc_unc, marker = 'D', markersize=7, linestyle='', capsize = 7, color='grey', label='avrg')

plt.xlabel('Beam energy (MeV)')
plt.ylabel('Beam current (nA)')
plt.legend()
plt.show()














#_____________Compartment 1_______________
# beam_current_list_comp1 = []
# beam_current_unc_list_comp1 = []

# beam_current_list_comp1.append(data_by_reaction['56CO']['beam_current'][0])
# beam_current_list_comp1.append(data_by_reaction['58CO']['beam_current'][0])
# beam_current_list_comp1.append(data_by_reaction['61CU']['beam_current'][0])
# beam_current_list_comp1.append(data_by_reaction['46SC']['beam_current'][0])
# beam_current_list_comp1.append(data_by_reaction['48V']['beam_current'][0])

# beam_current_unc_list_comp1.append(data_by_reaction['56CO']['beam_current_unc'][0])
# beam_current_unc_list_comp1.append(data_by_reaction['58CO']['beam_current_unc'][0])
# beam_current_unc_list_comp1.append(data_by_reaction['61CU']['beam_current_unc'][0])
# beam_current_unc_list_comp1.append(data_by_reaction['46SC']['beam_current_unc'][0])
# beam_current_unc_list_comp1.append(data_by_reaction['48V']['beam_current_unc'][0])

# avrg_bc_comp1, sigma_avrg_bc_comp1 = calculate_weighted_average_beam_current(beam_current_list_comp1, beam_current_unc_list_comp1)
# print(f'Avrg bc comp1: {avrg_bc_comp1} +- {sigma_avrg_bc_comp1}')


# #_____________Compartment 2_______________
# beam_current_list_comp2 = []
# beam_current_unc_list_comp2 = []

# beam_current_list_comp2.append(data_by_reaction['56CO']['beam_current'][1])
# beam_current_list_comp2.append(data_by_reaction['58CO']['beam_current'][1])
# beam_current_list_comp2.append(data_by_reaction['61CU']['beam_current'][1])
# beam_current_list_comp2.append(data_by_reaction['46SC']['beam_current'][1])
# beam_current_list_comp2.append(data_by_reaction['48V']['beam_current'][1])

# beam_current_unc_list_comp2.append(data_by_reaction['56CO']['beam_current_unc'][1])
# beam_current_unc_list_comp2.append(data_by_reaction['58CO']['beam_current_unc'][1])
# beam_current_unc_list_comp2.append(data_by_reaction['61CU']['beam_current_unc'][1])
# beam_current_unc_list_comp2.append(data_by_reaction['46SC']['beam_current_unc'][1])
# beam_current_unc_list_comp2.append(data_by_reaction['48V']['beam_current_unc'][1])

# avrg_bc_comp2, sigma_avrg_bc_comp2 = calculate_weighted_average_beam_current(beam_current_list_comp2, beam_current_unc_list_comp2)
# print(f'Avrg bc comp2: {avrg_bc_comp2} +- {sigma_avrg_bc_comp2}')


# #_____________Compartment 3_______________
# beam_current_list_comp3 = []
# beam_current_unc_list_comp3 = []

# beam_current_list_comp3.append(data_by_reaction['56CO']['beam_current'][2])
# beam_current_list_comp3.append(data_by_reaction['58CO']['beam_current'][2])
# beam_current_list_comp3.append(data_by_reaction['61CU']['beam_current'][2])
# beam_current_list_comp3.append(data_by_reaction['46SC']['beam_current'][2])
# beam_current_list_comp3.append(data_by_reaction['48V']['beam_current'][2])

# beam_current_unc_list_comp3.append(data_by_reaction['56CO']['beam_current_unc'][2])
# beam_current_unc_list_comp3.append(data_by_reaction['58CO']['beam_current_unc'][2])
# beam_current_unc_list_comp3.append(data_by_reaction['61CU']['beam_current_unc'][2])
# beam_current_unc_list_comp3.append(data_by_reaction['46SC']['beam_current_unc'][2])
# beam_current_unc_list_comp3.append(data_by_reaction['48V']['beam_current_unc'][2])

# avrg_bc_comp3, sigma_avrg_bc_comp3 = calculate_weighted_average_beam_current(beam_current_list_comp3, beam_current_unc_list_comp3)
# print(f'Avrg bc comp3: {avrg_bc_comp3} +- {sigma_avrg_bc_comp3}')


# #_____________Compartment 4_______________
# beam_current_list_comp4 = []
# beam_current_unc_list_comp4 = []

# beam_current_list_comp4.append(data_by_reaction['56CO']['beam_current'][3])
# beam_current_list_comp4.append(data_by_reaction['58CO']['beam_current'][3])
# beam_current_list_comp4.append(data_by_reaction['61CU']['beam_current'][3])
# beam_current_list_comp4.append(data_by_reaction['46SC']['beam_current'][3])
# beam_current_list_comp4.append(data_by_reaction['48V']['beam_current'][3])

# beam_current_unc_list_comp4.append(data_by_reaction['56CO']['beam_current_unc'][3])
# beam_current_unc_list_comp4.append(data_by_reaction['58CO']['beam_current_unc'][3])
# beam_current_unc_list_comp4.append(data_by_reaction['61CU']['beam_current_unc'][3])
# beam_current_unc_list_comp4.append(data_by_reaction['46SC']['beam_current_unc'][3])
# beam_current_unc_list_comp4.append(data_by_reaction['48V']['beam_current_unc'][3])

# avrg_bc_comp4, sigma_avrg_bc_comp4 = calculate_weighted_average_beam_current(beam_current_list_comp4, beam_current_unc_list_comp4)
# print(f'Avrg bc comp4: {avrg_bc_comp4} +- {sigma_avrg_bc_comp4}')


