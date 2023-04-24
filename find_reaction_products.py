
import numpy as np 


def calculate_Q_value(initial_mass_array, product_mass_array):
    #Function to calculate the Q-value of a reaction. Takes an array of the masses of all initial 
    #nuclei and an array of the masses of all the reaction products as input. Make sure that all 
    #the masses have the same unit (recommend MeV)
    Q_value = np.sum(initial_mass_array)-np.sum(product_mass_array)
    return Q_value



def find_possible_reactions(target_name, target_mass, Z_target, N_target, product_name, product_mass, Z_product, N_product, beamEnergy):

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
            print(f'{product_name} can be produced through the {target_name}(d,{-delta_N+1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,{-delta_N+1}n) . Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == 0 and delta_N <= 1): #(d,pXn)
        print('(d,pXn)')

        product_mass_array = np.array([product_mass, proton_mass, neutron_mass*(-delta_N+1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,p{-delta_N+1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,p{-delta_N+1}n). Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == 0 and delta_N <= 0): #(d,dXn)
        print('(d,dXn)')

        product_mass_array = np.array([product_mass, deuteron_mass, neutron_mass*(-delta_N)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,d{-delta_N}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,d{-delta_N}n). Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == 0 and delta_N <= -1): #(d,tXn)
        print('(d,tXn)')

        product_mass_array = np.array([product_mass, tritium_mass, neutron_mass*(-delta_N-1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,t{-delta_N-1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,t{-delta_N-1}n). Avaiable energy: {energyAvailable:.2f}MeV')


    if (delta_Z == -1 and delta_N <= -1): #(d,aXn)
        print('(d,aXn)')

        product_mass_array = np.array([product_mass, alpha_mass, neutron_mass*(-delta_N-1)])
        Q_value = calculate_Q_value(initial_mass_array, product_mass_array)
        energyAvailable = Q_value + beamEnergy

        if (energyAvailable>0):
            reaction_not_found = False
            print(f'{product_name} can be produced through the {target_name}(d,a{-delta_N-1}n){product_name} reaction.\nQ-value + beam energy = {energyAvailable:.2f}MeV\n')
        else:
            print(f'Not enough energy available to produce {product_name} through (d,a{-delta_N-1}n). Avaiable energy: {energyAvailable:.2f}MeV')

    if reaction_not_found:
        print('No possible reaction found')




if __name__ == '__main__':
    u = 931.49410242 #[MeV/c^2]

    Nb_isotopes = ['86Nb', '87Nb', '88Nb', '89Nb', '90Nb', '91Nb', '92Nb', '93Nb', '94Nb', '95Nb', '96Nb', '97Nb', '98Nb']
    Nb_masses = np.array([85.925781536, 86.920692473, 87.918226476, 88.913444696, 89.911259201, 90.906990256, 91.907188580, 92.906373170, 93.907279001, 94.906831110, 95.908101586, 96.908101622, 97.910332645]) #[u] from NNDC

    Zr_target_isotopes = ['90Zr', '91Zr', '92Zr', '93Zr', '94Zr', '96Zr']
    Zr_target_masses = np.array([89.904698755, 90.905640205, 91.905035336, 92.906470661, 93.906312523, 95.908277615]) #[u] from NNDC

    Zr_isotopes_for_d_pXn = ['86Zr', '87Zr', '88Zr', '89Zr', '90Zr', '91Zr', '92Zr', '93Zr', '94Zr', '95Zr', '96Zr', '97Zr']
    Zr_masses_for_d_pXn = np.array([85.916296814, 86.914817338, 87.910220715, 88.908879751, 89.904698755, 90.905640205, 91.905035336, 92.906470661, 93.906312523, 94.908040276, 95.908277615, 96.910963802]) #[u] from NNDC

    find_possible_reactions('90Zr', 89.904698755*u, 40, 50, '88Y', 87.909501274*u, 39, 49, 30)




    



# def Q_value_for_N_neutrons_in_d_Xn_reactions(Zr_mass, Nb_product_mass, N_neutrons):
#     #Function for calculating the Q-value for Zr(d,x*n)Nb for different Zr targets with different numbers of nutrons evaporating 
#     deuteron_mass = 1875.6 #[MeV/c^2]
#     neutron_mass = 939.5654133 #[MeV/c^2]

#     Q_value = Zr_mass + deuteron_mass - Nb_product_mass - N_neutrons*neutron_mass #[MeV]

#     return Q_value



# def Q_value_for_N_neutrons_in_d_pXn_reactions(Zr_mass, Zr_product_mass, N_neutrons):
#     deuteron_mass = 1875.6 #[MeV/c^2]
#     neutron_mass = 939.5654133 #[MeV/c^2]
#     proton_mass = 938.272088 #[MeV/c^2]

#     Q_value = Zr_mass + deuteron_mass - Zr_product_mass - proton_mass - N_neutrons*neutron_mass #[MeV]

#     return Q_value



# def Energy_available(Q_value, beamEnergy=30): 
#     energyAvailable = Q_value + beamEnergy

#     return energyAvailable



# def find_possible_reaction_products_from_d_Xn_reactions(Zr_isotopes, Zr_masses, Nb_isotopes, Nb_masses):

#     u = 931.49410242 #[MeV/c^2]
#     Zr_masses = Zr_masses*u #[MeV/c^2]
#     Nb_masses = Nb_masses*u #[MeV/c^2]

#     #Loop over every Zr isotopes with A=90-94:
#     for i in range(len(Zr_isotopes)-1): #For 90Zr-94Zr:

#         N_neutrons = 0
#         Q_value = Q_value_for_N_neutrons_in_d_Xn_reactions(np.abs(Zr_masses[i]), np.abs(Nb_masses[i+6-N_neutrons]), N_neutrons)
#         energyAvailable = Energy_available(Q_value)
#         print('\n')
#         print(f'Q-value for producing {Nb_isotopes[i+6-N_neutrons]} from {Zr_isotopes[i]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')

#         while energyAvailable > 0: #Increasing the number of neutrons evaporated until we get a negative Q-value

#             N_neutrons += 1
#             Q_value = Q_value_for_N_neutrons_in_d_Xn_reactions(np.abs(Zr_masses[i]), np.abs(Nb_masses[i+6-N_neutrons]), N_neutrons)
#             energyAvailable = Energy_available(Q_value)


#             if energyAvailable > 0:
#                 print(f'Q-value for producing {Nb_isotopes[i+6-N_neutrons]} from {Zr_isotopes[i]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')

#     # For 96Zr:
#     N_neutrons = 0
#     Q_value = Q_value_for_N_neutrons_in_d_Xn_reactions(np.abs(Zr_masses[-1]), np.abs(Nb_masses[len(Zr_isotopes)+6-N_neutrons]), N_neutrons)
#     energyAvailable = Energy_available(Q_value)

#     print('\n')
#     print(f'Q-value for producing {Nb_isotopes[len(Zr_isotopes)+6-N_neutrons]} from {Zr_isotopes[-1]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')


#     while energyAvailable > 0: #Increasing the number of neutrons evaporated until we get a negative Q-value

#         N_neutrons += 1
#         Q_value = Q_value_for_N_neutrons_in_d_Xn_reactions(np.abs(Zr_masses[-1]), np.abs(Nb_masses[len(Zr_isotopes)+6-N_neutrons]), N_neutrons)
#         energyAvailable = Energy_available(Q_value)


#         if energyAvailable > 0:
#             print(f'Q-value for producing {Nb_isotopes[len(Zr_isotopes)+6-N_neutrons]} from {Zr_isotopes[-1]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')




# def find_possible_reaction_products_from_d_pXn_reactions(Zr_target_isotopes, Zr_target_masses, Zr_prod_isotopes, Zr_prod_masses):

#     u = 931.49410242 #[MeV/c^2]
#     Zr_target_masses = Zr_target_masses*u #[MeV/c^2] 
#     Zr_prod_masses = Zr_prod_masses*u #[MeV/c^2]

#     #Loop over every Zr isotopes with A=90-94:
#     for i in range(len(Zr_target_isotopes)-1): #For 90Zr-94Zr:

#         N_neutrons = 0
#         Q_value = Q_value_for_N_neutrons_in_d_pXn_reactions(np.abs(Zr_target_masses[i]), np.abs(Zr_prod_masses[i+5-N_neutrons]), N_neutrons)
#         energyAvailable = Energy_available(Q_value)
#         print('\n')
#         print(f'Q-value for producing {Zr_prod_isotopes[i+5-N_neutrons]} from {Zr_target_isotopes[i]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')

#         while energyAvailable > 0: #Increasing the number of neutrons evaporated until we get a negative Q-value

#             N_neutrons += 1
#             Q_value = Q_value_for_N_neutrons_in_d_pXn_reactions(np.abs(Zr_target_masses[i]), np.abs(Zr_prod_masses[i+5-N_neutrons]), N_neutrons)
#             energyAvailable = Energy_available(Q_value)


#             if energyAvailable > 0:
#                 print(f'Q-value for producing {Zr_prod_isotopes[i+5-N_neutrons]} from {Zr_target_isotopes[i]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')

#     # # For 96Zr:
#     # N_neutrons = 0
#     # Q_value = Q_value_for_N_neutrons_in_d_pXn_reactions(np.abs(Zr_target_masses[-1]), np.abs(Zr_prod_masses[len(Zr_target_isotopes)-N_neutrons]), N_neutrons)
#     # energyAvailable = Energy_available(Q_value)

#     # print('\n')
#     # print(f'Q-value for producing {Zr_prod_isotopes[len(Zr_target_isotopes)-N_neutrons]} from {Zr_target_isotopes[-1]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')


#     # while energyAvailable > 0: #Increasing the number of neutrons evaporated until we get a negative Q-value

#     #     N_neutrons += 1
#     #     Q_value = Q_value_for_N_neutrons_in_d_pXn_reactions(np.abs(Zr_target_masses[-1]), np.abs(Zr_prod_masses[len(Zr_target_isotopes)-N_neutrons]), N_neutrons)
#     #     energyAvailable = Energy_available(Q_value)


#     #     if energyAvailable > 0:
#     #         print(f'Q-value for producing {Zr_prod_isotopes[len(Zr_target_isotopes)-N_neutrons]} from {Zr_target_isotopes[-1]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')




# if __name__ == '__main__':

#     #Mass excess for possible target and product nuclei 
#     Nb_isotopes = ['86Nb', '87Nb', '88Nb', '89Nb', '90Nb', '91Nb', '92Nb', '93Nb', '94Nb', '95Nb', '96Nb', '97Nb', '98Nb']
#     Nb_masses = np.array([85.925781536, 86.920692473, 87.918226476, 88.913444696, 89.911259201, 90.906990256, 91.907188580, 92.906373170, 93.907279001, 94.906831110, 95.908101586, 96.908101622, 97.910332645]) #[u] from NNDC

#     Zr_target_isotopes = ['90Zr', '91Zr', '92Zr', '93Zr', '94Zr', '96Zr']
#     Zr_target_masses = np.array([89.904698755, 90.905640205, 91.905035336, 92.906470661, 93.906312523, 95.908277615]) #[u] from NNDC

#     Zr_isotopes_for_d_pXn = ['86Zr', '87Zr', '88Zr', '89Zr', '90Zr', '91Zr', '92Zr', '93Zr', '94Zr', '95Zr', '96Zr', '97Zr']
#     Zr_masses_for_d_pXn = np.array([85.916296814, 86.914817338, 87.910220715, 88.908879751, 89.904698755, 90.905640205, 91.905035336, 92.906470661, 93.906312523, 94.908040276, 95.908277615, 96.910963802]) #[u] from NNDC


#     print('______________________________________________(d,Xn)__________________________________________________')
#     find_possible_reaction_products_from_d_Xn_reactions(Zr_target_isotopes, Zr_target_masses, Nb_isotopes, Nb_masses)
#     print('\n______________________________________________(d,pXn)__________________________________________________')
#     find_possible_reaction_products_from_d_pXn_reactions(Zr_target_isotopes, Zr_target_masses, Zr_isotopes_for_d_pXn, Zr_masses_for_d_pXn)

    
# if __name__ == '__main__':
#     #Mass excess for possible target and product nuclei 
#     Nb_isotopes = ['86Nb', '87Nb', '88Nb', '89Nb', '90Nb', '91Nb', '92Nb', '93Nb', '94Nb', '95Nb', '96Nb', '97Nb', '98Nb']
#     Nb_masses = np.array([85.925781536, 86.920692473, 87.918226476, 88.913444696, 89.911259201, 90.906990256, 91.907188580, 92.906373170, 93.907279001, 94.906831110, 95.908101586, 96.908101622, 97.910332645]) #[u] from NNDC

#     Zr_target_isotopes = ['90Zr', '91Zr', '92Zr', '93Zr', '94Zr', '96Zr']
#     Zr_target_masses = np.array([89.904698755, 90.905640205, 91.905035336, 92.906470661, 93.906312523, 95.908277615]) #[u] from NNDC

#     Zr_isotopes_for_d_pXn = ['86Zr', '87Zr', '88Zr', '89Zr', '90Zr', '91Zr', '92Zr', '93Zr', '94Zr', '95Zr', '96Zr', '97Zr']
#     Zr_masses_for_d_pXn = np.array([85.916296814, 86.914817338, 87.910220715, 88.908879751, 89.904698755, 90.905640205, 91.905035336, 92.906470661, 93.906312523, 94.908040276, 95.908277615, 96.910963802]) #[u] from NNDC


#     print('______________________________________________(d,Xn)__________________________________________________')
#     find_possible_reaction_products_from_d_Xn_reactions(Zr_target_isotopes, Zr_target_masses, Nb_isotopes, Nb_masses)
#     print('\n______________________________________________(d,pXn)__________________________________________________')
#     find_possible_reaction_products_from_d_pXn_reactions(Zr_target_isotopes, Zr_target_masses, Zr_isotopes_for_d_pXn, Zr_masses_for_d_pXn)




        

        
    

