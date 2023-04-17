
import numpy as np 



# def calc_Q_value(target_mass_excess, product_mass_excess, projectile_mass=1875.6,  beamEnergy=30):  #All input values must have units of MeV
#     #Calculating Q-value in MeV
#     energyStart = np.abs(target_mass_excess) + projectile_mass + beamEnergy
#     energyAfter = energyStart - np.abs(product_mass_excess)

#     return energyAfter


def Q_value_for_N_neutrons(Zr_mass, Nb_product_mass, N_neutrons, beamEnergy=30):
    deuteron_mass = 1875.6 #[MeV]
    neutron_mass = 939.5654133 #[MeV]

    Q_value = Zr_mass + deuteron_mass + beamEnergy - Nb_product_mass - N_neutrons*neutron_mass

    return Q_value



if __name__ == '__main__':

    Nb_isotopes = ['86Nb', '87Nb', '88Nb', '89Nb', '90Nb', '91Nb', '92Nb', '93Nb', '94Nb', '95Nb', '96Nb', '97Nb', '98Nb']
    Nb_mass_excess = np.array([-69134, -73874, -7.617E+4, -80626, -82662, -86638, -86453.3, -87212.8, -86369.1, -86786.3, -85602.83, -85603, -83525])/1000 #[MeV]

    Zr_isotopes = ['90Zr', '91Zr', '92Zr','93Zr', '94Zr', '96Zr']
    Zr_mass_excess = np.array([-88772.55, -87895.59, -88459.02, -87122.0, -87269.33, -85438.86])/1000 #[MeV]


    #for Zr = i har vi Nb-compound = i+6
    # Nb product vil være Nb-compund - N_neutrons = i+6-N_neutrons

    for i in range(len(Zr_isotopes)-1): #For 90Zr-94Zr:

        N_neutrons = 0
        Q_value = Q_value_for_N_neutrons(np.abs(Zr_mass_excess[i]), np.abs(Nb_mass_excess[i+6-N_neutrons]), N_neutrons)
        print('\n')
        print(f'Q-value for producing {Nb_isotopes[i+6-N_neutrons]} form {Zr_isotopes[i]} is {Q_value:.2f}MeV.')

        while Q_value > 0:

            N_neutrons += 1
            Q_value = Q_value_for_N_neutrons(np.abs(Zr_mass_excess[i]), np.abs(Nb_mass_excess[i+6-N_neutrons]), N_neutrons)

            if Q_value > 0:
                print(f'Q-value for producing {Nb_isotopes[i+6-N_neutrons]} form {Zr_isotopes[i]} is {Q_value:.2f}MeV.')

    # For 96Zr:
    N_neutrons = 0
    Q_value = Q_value_for_N_neutrons(np.abs(Zr_mass_excess[-1]), np.abs(Nb_mass_excess[len(Zr_isotopes)+6-N_neutrons]), N_neutrons)
    print('\n')
    print(f'Q-value for producing {Nb_isotopes[len(Zr_isotopes)+6-N_neutrons]} form {Zr_isotopes[-1]} is {Q_value:.2f}MeV.')


    while Q_value > 0:

        N_neutrons += 1
        Q_value = Q_value_for_N_neutrons(np.abs(Zr_mass_excess[-1]), np.abs(Nb_mass_excess[len(Zr_isotopes)+6-N_neutrons]), N_neutrons)

        if Q_value > 0:
            print(f'Q-value for producing {Nb_isotopes[len(Zr_isotopes)+6-N_neutrons]} form {Zr_isotopes[-1]} is {Q_value:.2f}MeV.')





        

        
    

