
import numpy as np 



def Q_value_for_N_neutrons(Zr_mass, Nb_product_mass, N_neutrons, beamEnergy=30):
    #Function for calculating the Q-value for Zr(d,x*n)Nb for different Zr targets with different numbers of nutrons evaporating 
    deuteron_mass = 1875.6 #[MeV/c^2]
    neutron_mass = 939.5654133 #[MeV/c^2]

    Q_value = Zr_mass + deuteron_mass - Nb_product_mass - N_neutrons*neutron_mass #[MeV]

    return Q_value



def Energy_available(Q_value, beamEnergy=30): 
    energyAvailable = Q_value + beamEnergy

    return energyAvailable



def find_possible_reaction_products_from_Zr_d_n_Nb_reactions(Zr_isotopes, Zr_masses, Nb_isotopes, Nb_masses):

    u = 931.49410242 #[MeV/c^2]
    Zr_masses = Zr_masses*u
    Nb_masses = Nb_masses*u

    #Loop over every Zr isotopes with A=90-94:
    for i in range(len(Zr_isotopes)-1): #For 90Zr-94Zr:

        N_neutrons = 0
        Q_value = Q_value_for_N_neutrons(np.abs(Zr_masses[i]), np.abs(Nb_masses[i+6-N_neutrons]), N_neutrons)
        energyAvailable = Energy_available(Q_value)
        print('\n')
        print(f'Q-value for producing {Nb_isotopes[i+6-N_neutrons]} from {Zr_isotopes[i]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')

        while energyAvailable > 0: #Increasing the number of neutrons evaporated until we get a negative Q-value

            N_neutrons += 1
            Q_value = Q_value_for_N_neutrons(np.abs(Zr_masses[i]), np.abs(Nb_masses[i+6-N_neutrons]), N_neutrons)
            energyAvailable = Energy_available(Q_value)


            if energyAvailable > 0:
                print(f'Q-value for producing {Nb_isotopes[i+6-N_neutrons]} from {Zr_isotopes[i]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')

    # For 96Zr:
    N_neutrons = 0
    Q_value = Q_value_for_N_neutrons(np.abs(Zr_masses[-1]), np.abs(Nb_masses[len(Zr_isotopes)+6-N_neutrons]), N_neutrons)
    energyAvailable = Energy_available(Q_value)

    print('\n')
    print(f'Q-value for producing {Nb_isotopes[len(Zr_isotopes)+6-N_neutrons]} from {Zr_isotopes[-1]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')


    while energyAvailable > 0: #Increasing the number of neutrons evaporated until we get a negative Q-value

        N_neutrons += 1
        Q_value = Q_value_for_N_neutrons(np.abs(Zr_masses[-1]), np.abs(Nb_masses[len(Zr_isotopes)+6-N_neutrons]), N_neutrons)
        energyAvailable = Energy_available(Q_value)


        if energyAvailable > 0:
            print(f'Q-value for producing {Nb_isotopes[len(Zr_isotopes)+6-N_neutrons]} from {Zr_isotopes[-1]} is {Q_value:.2f}MeV, and the available energy after the reaction is {energyAvailable:.2f}MeV.')




if __name__ == '__main__':

    #Mass excess for possible target and product nuclei 
    Nb_isotopes = ['86Nb', '87Nb', '88Nb', '89Nb', '90Nb', '91Nb', '92Nb', '93Nb', '94Nb', '95Nb', '96Nb', '97Nb', '98Nb']
    Nb_masses = np.array([85.925781536, 86.920692473, 87.918226476, 88.913444696, 89.911259201, 90.906990256, 91.907188580, 92.906373170, 93.907279001, 94.906831110, 95.908101586, 96.908101622, 97.910332645]) #[u] from NNDC

    Zr_isotopes = ['90Zr', '91Zr', '92Zr', '93Zr', '94Zr', '96Zr']
    Zr_masses = np.array([89.904698755, 90.905640205, 91.905035336, 92.906470661, 93.906312523, 95.908277615]) #[u] from NNDC


    find_possible_reaction_products_from_Zr_d_n_Nb_reactions(Zr_isotopes, Zr_masses, Nb_isotopes, Nb_masses)

    





        

        
    

