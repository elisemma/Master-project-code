import numpy as np 


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






data_by_reaction = {'56CO': {'beam_current': [91.11279483618009, 95.48943423201047, 96.67759656591667, 87.90385977930694, 79.76353046504703], 
                             'beam_current_unc': [17.764379182058054, 17.550901526554565, 16.71724750021088, 13.273019225954833, 9.946307812204912], 
                             'energy': [48.25110370000002, 41.876494533333336, 37.24984646666667, 32.115906233333334, 26.2370075]}, 
                    '46SC': {'beam_current': [98.38985492917791, 92.43719344140617, 94.87157275894245, 92.64152715822723, 89.76044811394085], 
                             'beam_current_unc': [9.06208185316921, 8.211563169758909, 8.188146100794526, 8.695843008929224, 7.418752250751914], 
                             'energy': [47.13022093333334, 40.62009246666667, 35.86358033333333, 30.5522245, 24.390982466666664]}, 
                    '48V':  {'beam_current': [88.69779045350535, 89.04935458091553, 92.2968314456004, 90.59958822021765, 89.8492617019322], 
                             'beam_current_unc': [12.20254718194584, 9.989933425556139, 9.93150127659525, 10.305181149242912, 9.318875229511063], 
                             'energy': [47.13022093333334, 40.62009246666667, 35.86358033333333, 30.5522245, 24.390982466666664]}}





#_____________Compartment 1_______________
beam_current_list_comp1 = []
beam_current_unc_list_comp1 = []

beam_current_list_comp1.append(data_by_reaction['56CO']['beam_current'][0])
beam_current_list_comp1.append(data_by_reaction['46SC']['beam_current'][0])
beam_current_list_comp1.append(data_by_reaction['48V']['beam_current'][0])

beam_current_unc_list_comp1.append(data_by_reaction['56CO']['beam_current_unc'][0])
beam_current_unc_list_comp1.append(data_by_reaction['46SC']['beam_current_unc'][0])
beam_current_unc_list_comp1.append(data_by_reaction['48V']['beam_current_unc'][0])

avrg_bc_comp1, sigma_avrg_bc_comp1 = calculate_weighted_average_beam_current(beam_current_list_comp1, beam_current_unc_list_comp1)
print(f'Avrg bc comp1: {avrg_bc_comp1} +- {sigma_avrg_bc_comp1}')


#_____________Compartment 2_______________
beam_current_list_comp2 = []
beam_current_unc_list_comp2 = []

beam_current_list_comp2.append(data_by_reaction['56CO']['beam_current'][1])
beam_current_list_comp2.append(data_by_reaction['46SC']['beam_current'][1])
beam_current_list_comp2.append(data_by_reaction['48V']['beam_current'][1])

beam_current_unc_list_comp2.append(data_by_reaction['56CO']['beam_current_unc'][1])
beam_current_unc_list_comp2.append(data_by_reaction['46SC']['beam_current_unc'][1])
beam_current_unc_list_comp2.append(data_by_reaction['48V']['beam_current_unc'][1])

avrg_bc_comp2, sigma_avrg_bc_comp2 = calculate_weighted_average_beam_current(beam_current_list_comp2, beam_current_unc_list_comp2)
print(f'Avrg bc comp2: {avrg_bc_comp2} +- {sigma_avrg_bc_comp2}')


#_____________Compartment 3_______________
beam_current_list_comp3 = []
beam_current_unc_list_comp3 = []

beam_current_list_comp3.append(data_by_reaction['56CO']['beam_current'][2])
beam_current_list_comp3.append(data_by_reaction['46SC']['beam_current'][2])
beam_current_list_comp3.append(data_by_reaction['48V']['beam_current'][2])

beam_current_unc_list_comp3.append(data_by_reaction['56CO']['beam_current_unc'][2])
beam_current_unc_list_comp3.append(data_by_reaction['46SC']['beam_current_unc'][2])
beam_current_unc_list_comp3.append(data_by_reaction['48V']['beam_current_unc'][2])

avrg_bc_comp3, sigma_avrg_bc_comp3 = calculate_weighted_average_beam_current(beam_current_list_comp3, beam_current_unc_list_comp3)
print(f'Avrg bc comp3: {avrg_bc_comp3} +- {sigma_avrg_bc_comp3}')


#_____________Compartment 4_______________
beam_current_list_comp4 = []
beam_current_unc_list_comp4 = []

beam_current_list_comp4.append(data_by_reaction['56CO']['beam_current'][3])
beam_current_list_comp4.append(data_by_reaction['46SC']['beam_current'][3])
beam_current_list_comp4.append(data_by_reaction['48V']['beam_current'][3])

beam_current_unc_list_comp4.append(data_by_reaction['56CO']['beam_current_unc'][3])
beam_current_unc_list_comp4.append(data_by_reaction['46SC']['beam_current_unc'][3])
beam_current_unc_list_comp4.append(data_by_reaction['48V']['beam_current_unc'][3])

avrg_bc_comp4, sigma_avrg_bc_comp4 = calculate_weighted_average_beam_current(beam_current_list_comp4, beam_current_unc_list_comp4)
print(f'Avrg bc comp4: {avrg_bc_comp4} +- {sigma_avrg_bc_comp4}')


#_____________Compartment 5_______________
beam_current_list_comp5 = []
beam_current_unc_list_comp5 = []

beam_current_list_comp5.append(data_by_reaction['56CO']['beam_current'][4])
beam_current_list_comp5.append(data_by_reaction['46SC']['beam_current'][4])
beam_current_list_comp5.append(data_by_reaction['48V']['beam_current'][4])

beam_current_unc_list_comp5.append(data_by_reaction['56CO']['beam_current_unc'][4])
beam_current_unc_list_comp5.append(data_by_reaction['46SC']['beam_current_unc'][4])
beam_current_unc_list_comp5.append(data_by_reaction['48V']['beam_current_unc'][4])

avrg_bc_comp5, sigma_avrg_bc_comp5 = calculate_weighted_average_beam_current(beam_current_list_comp5, beam_current_unc_list_comp5)
print(f'Avrg bc comp5: {avrg_bc_comp5} +- {sigma_avrg_bc_comp5}')

