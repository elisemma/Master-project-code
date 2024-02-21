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





#_____________Compartment 1_______________
beam_current_list_comp1 = []
beam_current_unc_list_comp1 = []

beam_current_list_comp1.append(data_by_reaction['56CO']['beam_current'][0])
beam_current_list_comp1.append(data_by_reaction['58CO']['beam_current'][0])
beam_current_list_comp1.append(data_by_reaction['61CU']['beam_current'][0])
beam_current_list_comp1.append(data_by_reaction['46SC']['beam_current'][0])
beam_current_list_comp1.append(data_by_reaction['48V']['beam_current'][0])

beam_current_unc_list_comp1.append(data_by_reaction['56CO']['beam_current_unc'][0])
beam_current_unc_list_comp1.append(data_by_reaction['58CO']['beam_current_unc'][0])
beam_current_unc_list_comp1.append(data_by_reaction['61CU']['beam_current_unc'][0])
beam_current_unc_list_comp1.append(data_by_reaction['46SC']['beam_current_unc'][0])
beam_current_unc_list_comp1.append(data_by_reaction['48V']['beam_current_unc'][0])

avrg_bc_comp1, sigma_avrg_bc_comp1 = calculate_weighted_average_beam_current(beam_current_list_comp1, beam_current_unc_list_comp1)
print(f'Avrg bc comp1: {avrg_bc_comp1} +- {sigma_avrg_bc_comp1}')


#_____________Compartment 2_______________
beam_current_list_comp2 = []
beam_current_unc_list_comp2 = []

beam_current_list_comp2.append(data_by_reaction['56CO']['beam_current'][1])
beam_current_list_comp2.append(data_by_reaction['58CO']['beam_current'][1])
beam_current_list_comp2.append(data_by_reaction['61CU']['beam_current'][1])
beam_current_list_comp2.append(data_by_reaction['46SC']['beam_current'][1])
beam_current_list_comp2.append(data_by_reaction['48V']['beam_current'][1])

beam_current_unc_list_comp2.append(data_by_reaction['56CO']['beam_current_unc'][1])
beam_current_unc_list_comp2.append(data_by_reaction['58CO']['beam_current_unc'][1])
beam_current_unc_list_comp2.append(data_by_reaction['61CU']['beam_current_unc'][1])
beam_current_unc_list_comp2.append(data_by_reaction['46SC']['beam_current_unc'][1])
beam_current_unc_list_comp2.append(data_by_reaction['48V']['beam_current_unc'][1])

avrg_bc_comp2, sigma_avrg_bc_comp2 = calculate_weighted_average_beam_current(beam_current_list_comp2, beam_current_unc_list_comp2)
print(f'Avrg bc comp2: {avrg_bc_comp2} +- {sigma_avrg_bc_comp2}')


#_____________Compartment 3_______________
beam_current_list_comp3 = []
beam_current_unc_list_comp3 = []

beam_current_list_comp3.append(data_by_reaction['56CO']['beam_current'][2])
beam_current_list_comp3.append(data_by_reaction['58CO']['beam_current'][2])
beam_current_list_comp3.append(data_by_reaction['61CU']['beam_current'][2])
beam_current_list_comp3.append(data_by_reaction['46SC']['beam_current'][2])
beam_current_list_comp3.append(data_by_reaction['48V']['beam_current'][2])

beam_current_unc_list_comp3.append(data_by_reaction['56CO']['beam_current_unc'][2])
beam_current_unc_list_comp3.append(data_by_reaction['58CO']['beam_current_unc'][2])
beam_current_unc_list_comp3.append(data_by_reaction['61CU']['beam_current_unc'][2])
beam_current_unc_list_comp3.append(data_by_reaction['46SC']['beam_current_unc'][2])
beam_current_unc_list_comp3.append(data_by_reaction['48V']['beam_current_unc'][2])

avrg_bc_comp3, sigma_avrg_bc_comp3 = calculate_weighted_average_beam_current(beam_current_list_comp3, beam_current_unc_list_comp3)
print(f'Avrg bc comp3: {avrg_bc_comp3} +- {sigma_avrg_bc_comp3}')


#_____________Compartment 4_______________
beam_current_list_comp4 = []
beam_current_unc_list_comp4 = []

beam_current_list_comp4.append(data_by_reaction['56CO']['beam_current'][3])
beam_current_list_comp4.append(data_by_reaction['58CO']['beam_current'][3])
beam_current_list_comp4.append(data_by_reaction['61CU']['beam_current'][3])
beam_current_list_comp4.append(data_by_reaction['46SC']['beam_current'][3])
beam_current_list_comp4.append(data_by_reaction['48V']['beam_current'][3])

beam_current_unc_list_comp4.append(data_by_reaction['56CO']['beam_current_unc'][3])
beam_current_unc_list_comp4.append(data_by_reaction['58CO']['beam_current_unc'][3])
beam_current_unc_list_comp4.append(data_by_reaction['61CU']['beam_current_unc'][3])
beam_current_unc_list_comp4.append(data_by_reaction['46SC']['beam_current_unc'][3])
beam_current_unc_list_comp4.append(data_by_reaction['48V']['beam_current_unc'][3])

avrg_bc_comp4, sigma_avrg_bc_comp4 = calculate_weighted_average_beam_current(beam_current_list_comp4, beam_current_unc_list_comp4)
print(f'Avrg bc comp4: {avrg_bc_comp4} +- {sigma_avrg_bc_comp4}')


