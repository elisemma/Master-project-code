
import numpy as np 


def calculate_Q_value(initial_mass_array, product_mass_array):
    #Function to calculate the Q-value of a reaction. Takes an array of the masses of all initial 
    #nuclei and an array of the masses of all the reaction products as input. Make sure that all 
    #the masses have the same unit (recommend MeV)
    Q_value = np.sum(initial_mass_array)-np.sum(product_mass_array)
    return Q_value



def find_possible_reactions_verbose(target_name, target_mass, Z_target, N_target, product_name, product_mass, Z_product, N_product, beamEnergy):
    #This function prints put a lot of output and is therefore nice to use if you want to take a deeper look into the channels for producing one specific nuclei

    delta_Z = Z_product-Z_target
    delta_N = N_product-N_target

    neutron_mass = 939.5654133 #[MeV/c^2]
    proton_mass = 938.272088 #[MeV/c^2]
    deuteron_mass = 1875.6 #[MeV/c^2]
    tritium_mass = 2809.432117 #[MeV/c^2]
    alpha_mass = 3727.3794066 #[MeV/c^2]

    initial_mass_array = np.array([target_mass, deuteron_mass])

    reaction_not_found = True

    if (delta_Z == 1 and delta_N <= 1): #(d,Xn)
        print('(d,Xn)')

        product_mass_array = np.array([product_mass, neutron_mass*(-delta_N+1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,{-delta_N+1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\nQ-value = {Q_value:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,{-delta_N+1}n) . Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == 0 and delta_N <= 1): #(d,pXn)
        print('(d,pXn)')

        product_mass_array = np.array([product_mass, proton_mass, neutron_mass*(-delta_N+1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,p{-delta_N+1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\nQ-value = {Q_value:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,p{-delta_N+1}n). Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == 0 and delta_N <= 0): #(d,dXn)
        print('(d,dXn)')

        product_mass_array = np.array([product_mass, deuteron_mass, neutron_mass*(-delta_N)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,d{-delta_N}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\nQ-value = {Q_value:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,d{-delta_N}n). Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == 0 and delta_N <= -1): #(d,tXn)
        print('(d,tXn)')

        product_mass_array = np.array([product_mass, tritium_mass, neutron_mass*(-delta_N-1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,t{-delta_N-1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\nQ-value = {Q_value:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,t{-delta_N-1}n). Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == -1 and delta_N <= -1): #(d,aXn)
        print('(d,aXn)')

        product_mass_array = np.array([product_mass, alpha_mass, neutron_mass*(-delta_N-1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,a{-delta_N-1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\nQ-value = {Q_value:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,a{-delta_N-1}n). Avaiable energy: {energyAvailable:.2f}MeV')

    if reaction_not_found:
        print('No possible reaction found')




def find_possible_reactions(target_name, target_mass, Z_target, N_target, product_name, product_mass, Z_product, N_product, beamEnergy):
    #This function does the same as find_possible_reactions_verbose, but does not print anything and returns a list of the possible reactions instead. This is nice when looping over many nuclei

    delta_Z = Z_product-Z_target
    delta_N = N_product-N_target

    neutron_mass = 939.5654133 #[MeV/c^2]
    proton_mass = 938.272088 #[MeV/c^2]
    deuteron_mass = 1875.6 #[MeV/c^2]
    tritium_mass = 2809.432117 #[MeV/c^2]
    alpha_mass = 3727.3794066 #[MeV/c^2]

    initial_mass_array = np.array([target_mass, deuteron_mass])

    ListOfPossibleReactions = []


    if (delta_Z == 1 and delta_N <= 1): #(d,Xn)
    
        product_mass_array = np.array([product_mass, neutron_mass*(-delta_N+1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            ListOfPossibleReactions.append(f'{target_name}(d,{-delta_N+1}n){product_name}')


    if (delta_Z == 0 and delta_N <= 1): #(d,pXn)

        product_mass_array = np.array([product_mass, proton_mass, neutron_mass*(-delta_N+1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            ListOfPossibleReactions.append(f'{target_name}(d,p{-delta_N+1}n){product_name}')
          

    if (delta_Z == 0 and delta_N <= 0): #(d,dXn)
        
        product_mass_array = np.array([product_mass, deuteron_mass, neutron_mass*(-delta_N)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            ListOfPossibleReactions.append(f'{target_name}(d,d{-delta_N}n){product_name}')
        

    if (delta_Z == 0 and delta_N <= -1): #(d,tXn)
        
        product_mass_array = np.array([product_mass, tritium_mass, neutron_mass*(-delta_N-1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            ListOfPossibleReactions.append(f'{target_name}(d,t{-delta_N-1}n){product_name}')
         

    if (delta_Z == -1 and delta_N <= -1): #(d,aXn)
        
        product_mass_array = np.array([product_mass, alpha_mass, neutron_mass*(-delta_N-1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            ListOfPossibleReactions.append(f'{target_name}(d,a{-delta_N-1}n){product_name}')
            

    return ListOfPossibleReactions




def iterating_over_targets_and_products(target_dict, poroduct_dict):
    #lopps over many nuclei to see which nuclei can be produced and prints out the possible reaction to produce each nuclei

    for product_name, product_info in poroduct_dict.items():
        ListOfListsOfPossibleReactions = []
       
        for target_name, target_info in target_dict.items():
            ListOfPossibleReactions = find_possible_reactions(target_name, target_info['mass']*u, target_info['Z'], target_info['N'], product_name, product_info['mass']*u, product_info['Z'], product_info['N'], 30)
            ListOfListsOfPossibleReactions.append(ListOfPossibleReactions)

        print(f'\nPossible reactions for producing {product_name}:')

        for list in ListOfListsOfPossibleReactions:
            for reaction in list:
                print(reaction)





if __name__ == '__main__':


    u = 931.49410242 #[MeV/c^2]
    target_dict = {'90Zr': {'mass': 89.904698755, 'Z': 40, 'N': 50}, 
                   '91Zr': {'mass': 90.905640205, 'Z': 40, 'N': 51}, 
                   '92Zr': {'mass': 91.905035336, 'Z': 40, 'N': 52}, 
                   '94Zr': {'mass': 93.906312523, 'Z': 40, 'N': 54}, 
                   '96Zr': {'mass': 95.908277615, 'Z': 40, 'N': 56}} #masses in [u] from NNDC

    product_dict = {'86Nb': {'mass': 85.925781536, 'Z': 41, 'N': 45},
                     '87Nb': {'mass': 86.920692473, 'Z': 41, 'N': 46},
                     '88Nb': {'mass': 87.918226476, 'Z': 41, 'N': 47},
                     '89Nb': {'mass': 88.913444696, 'Z': 41, 'N': 48},
                     '90Nb': {'mass': 89.911259201, 'Z': 41, 'N': 49},
                     '91Nb': {'mass': 90.906990256, 'Z': 41, 'N': 50},
                     '92Nb': {'mass': 91.907188580, 'Z': 41, 'N': 51},
                     '93Nb': {'mass': 92.906373170, 'Z': 41, 'N': 52},
                     '94Nb': {'mass': 93.907279001, 'Z': 41, 'N': 53},
                     '95Nb': {'mass': 94.906831110, 'Z': 41, 'N': 54},
                     '96Nb': {'mass': 95.908101586, 'Z': 41, 'N': 55},
                     '97Nb': {'mass': 96.908101622, 'Z': 41, 'N': 56},
                     '98Nb': {'mass': 97.910332645, 'Z': 41, 'N': 57},
 
                     '90Zr': {'mass': 89.904698755, 'Z': 40, 'N': 50}, 
                     '91Zr': {'mass': 90.905640205, 'Z': 40, 'N': 51}, 
                     '92Zr': {'mass': 91.905035336, 'Z': 40, 'N': 52}, 
                     '93Zr': {'mass': 92.906470661, 'Z': 40, 'N': 53},
                     '94Zr': {'mass': 93.906312523, 'Z': 40, 'N': 54}, 
                     '95Zr': {'mass': 94.908040276, 'Z': 40, 'N': 55},
                     '96Zr': {'mass': 95.908277615, 'Z': 40, 'N': 56},

                     '86Y': {'mass': 85.914886095, 'Z': 39, 'N': 47},
                     '87Y': {'mass': 86.910876100, 'Z': 39, 'N': 48},
                     '88Y': {'mass': 87.909501274, 'Z': 39, 'N': 49},
                     '89Y': {'mass': 88.905838156, 'Z': 39, 'N': 50},
                     '90Y': {'mass': 89.907141749, 'Z': 39, 'N': 51},
                     '91Y': {'mass': 90.907298048, 'Z': 39, 'N': 52},
                     '92Y': {'mass': 91.908945752, 'Z': 39, 'N': 53},
                     '93Y': {'mass': 92.909578434, 'Z': 39, 'N': 54},
                     '94Y': {'mass': 93.911592062, 'Z': 39, 'N': 55}} #masses in [u] from NNDC

    

iterating_over_targets_and_products(target_dict, product_dict)

print('\n')

find_possible_reactions_verbose('96Zr', target_dict['96Zr']['mass']*u, target_dict['96Zr']['Z'], target_dict['96Zr']['N'], 
                                '94Y', product_dict['94Y']['mass']*u, product_dict['94Y']['Z'], product_dict['94Y']['N'], 30)
    

