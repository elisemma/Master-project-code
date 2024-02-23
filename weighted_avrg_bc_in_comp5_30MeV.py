import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge

#________________________________________________________________________30 MeV ______________________________________________________________________________________________
number_of_monitor_foils_30MeV = 1
monitor_reactions_per_foil_30MeV = np.array([2]) 
number_of_comp_30MeV = 1
number_of_monitor_reactions_30MeV = 2

A0_30MeV = np.array([ # [Bq]
    [25.452633286057697*12.642, 461.7853246484494*12.642]     # Comp 5 activities for 58Co, 61Cu
])  
  
sigma_A0_30MeV = np.array([ # [Bq]
    [2.465479823966241*12.642,  118.4114034137288*12.642]     # Comp 5 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
])  


mass_density_30MeV = np.array([ # [nuclei/cm^2]                    
    [2.333816317282693e+20,  2.333816317282693e+20]   # Comp 5 mass dens for Ni05, Ni05, Ni05, Ti05, Ti05
])  
  
sigma_mass_density_30MeV = np.array([ # [nuclei/cm^2]  
    [1.9488015308666957e+17, 1.9488015308666957e+17]   # Comp 5 mass dens unc for Ni05, Ni05, Ni05, Ti05, Ti05
])   
  
lambda__30MeV = np.log(2)/np.array([ # [1/s]                      
    [70.86*24*3600, 3.339*3600]     # Comp 5 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
reaction_integral_30MeV = 1e4*np.array([ # [cm^2]
    [3.3735628091502424e-31, 8.16618313949607e-31]   # Comp 5 int for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
uncertainty_integral_30MeV = 1e4*np.array([ # [cm^2]
    [3.916078608773901e-32,  9.173104893466661e-32]  # Comp 5 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
irr_time_30MeV = np.ones((number_of_comp_30MeV, number_of_monitor_reactions_30MeV))*1200 # [s]
sigma_irr_time_30MeV = np.ones((number_of_comp_30MeV, number_of_monitor_reactions_30MeV))*3 # [s]




A0 = A0_30MeV
sigma_A0 = sigma_A0_30MeV 
mass_density = mass_density_30MeV
sigma_mass_density = sigma_mass_density_30MeV 
lambda_= lambda__30MeV 
reaction_integral = reaction_integral_30MeV
uncertainty_integral = uncertainty_integral_30MeV 
irr_time = irr_time_30MeV 
sigma_irr_time = sigma_irr_time_30MeV
number_of_monitor_foils = number_of_monitor_foils_30MeV
monitor_reactions_per_foil = monitor_reactions_per_foil_30MeV











def Average_BeamCurrent(A0, sigma_A0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time, csv_filename='averaged_currents.csv', save_csv=False):


    # def decomment(csvfile):
    #     for row in csvfile:
    #         raw = row.split('#')[0].strip()
    #         if raw: yield raw

    # def read_csv(name_of_csv_file):
    #     results = []
    #     with open(name_of_csv_file) as csvfile:
    #         reader = csv.reader(decomment(csvfile))
    #         for row in reader:
    #             results.append(row)

    #     return np.asarray(results, dtype=float)


    def beam_current(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        return (elementary_charge*1e9 * A0) / (rho_dr * (1 - np.exp(-lambdas*t_irradiation)) * reaction_integral)


    # Numerical partial derivatives
    def dIdA0(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * A0
        return ((beam_current(A0 + (delta_x/2), rho_dr, lambdas, t_irradiation, reaction_integral) - beam_current(A0 - (delta_x/2), rho_dr, lambdas, t_irradiation, reaction_integral)) / delta_x)
    def dIdRhoDr(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * rho_dr
        return ((beam_current(A0, rho_dr + (delta_x/2), lambdas, t_irradiation, reaction_integral) - beam_current(A0, rho_dr - (delta_x/2), lambdas, t_irradiation, reaction_integral)) / delta_x)
    def dIdLambda(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * lambdas
        return ((beam_current(A0, rho_dr, lambdas + (delta_x/2), t_irradiation, reaction_integral) - beam_current(A0, rho_dr, lambdas - (delta_x/2), t_irradiation, reaction_integral)) / delta_x)
    def dIdTIrradiation(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * t_irradiation
        return ((beam_current(A0, rho_dr, lambdas, t_irradiation + (delta_x/2), reaction_integral) - beam_current(A0, rho_dr, lambdas, t_irradiation - (delta_x/2), reaction_integral)) / delta_x)
    def dIdIntegral(A0, rho_dr, lambdas, t_irradiation, reaction_integral):
        delta_x = 1E-8 * reaction_integral
        return ((beam_current(A0, rho_dr, lambdas, t_irradiation, reaction_integral + (delta_x/2)) - beam_current(A0, rho_dr, lambdas, t_irradiation, reaction_integral - (delta_x/2))) / delta_x)


    # Approximate uncertainties in beam current
    def sigma_I_approximate(A0, rho_dr, lambdas, t_irradiation, reaction_integral, unc_A0, unc_rho_dr, unc_lambdas, unc_t_irradiation, unc_reaction_integral):
        approx_error = np.sqrt(np.power(dIdA0(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_A0,2) +
        np.power(dIdRhoDr(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_rho_dr,2) +
        np.power(dIdLambda(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_lambdas,2) +
        np.power(dIdTIrradiation(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_t_irradiation,2) +
        #np.power( unc_reaction_integral,2))
        np.power(dIdIntegral(A0, rho_dr, lambdas, t_irradiation, reaction_integral) * unc_reaction_integral,2))
        # approx_error = 0
        # approx_error = np.power(dIdA0(A0, rho_dr, lambdas, t_irradiation, reaction_integral),2)
        return approx_error


    number_of_monitor_reactions = 0
    #print('yo')
    submatrix_lower_indices = np.zeros(number_of_monitor_foils)
    submatrix_upper_indices = np.zeros(number_of_monitor_foils)


    #### PARAMETERS IN THE FUNCTION des19_BeamCurrent.py
    #A0, dA0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time



    # Load in activation data
    #
    # All in nuclei / cm^2
    #                       Cu       Ti
    areal_density = mass_density
    uncertainty_areal_density = sigma_mass_density

    # All in Bq
    #                      Sc46    V48    Zn62    Zn63
    EoB_activities = A0
    uncertainty_EoB_activities = sigma_A0
    # Normalized integral(sigma * dPhidE * dE)
    #
    reaction_integral = reaction_integral
    #unc_rxn_integral = uncertainty_integral #rxn_int *     percent_rn_uncertainties      #rxn = reactions
    unc_rxn_integral = uncertainty_integral * reaction_integral#rxn_int *   percent_rn_uncertainties      #rxn = reactions
    # All in 1/s
    lambdas = lambda_
    # print(type(lambdas))
    # All in s
    t_irradiation = irr_time
    uncertainty_t_irradiation = sigma_irr_time


    for i in range(0, number_of_monitor_foils):
        submatrix_lower_indices[i] = number_of_monitor_reactions
        number_of_monitor_reactions += monitor_reactions_per_foil[i]
        submatrix_upper_indices[i] = number_of_monitor_reactions

    submatrix_lower_indices=submatrix_lower_indices.astype(int)
    submatrix_upper_indices=submatrix_upper_indices.astype(int)


    # Set up correlation matrices
    # Lambda is completely uncorrelated, except on the diagonal
    corr_lambda = np.zeros((number_of_monitor_reactions,number_of_monitor_reactions))
    # Areal density is completely uncorrelated, except within one foil's submatrix (fully correlated there)
    corr_areal_density = np.zeros((number_of_monitor_reactions,number_of_monitor_reactions))
    # Reaction integral is completely uncorrelated, except within one foil's submatrix (30% correlated there)
    corr_reaction_integral = np.zeros((number_of_monitor_reactions,number_of_monitor_reactions))
    # Irradiation length is completely correlated (same for all foils)
    corr_t_irradiation = np.ones((number_of_monitor_reactions,number_of_monitor_reactions))
    # EoB activitoes are partially uncorrelated (similar subset of efficiencies)
    corr_EoB_activities = 0.3 * np.ones((number_of_monitor_reactions,number_of_monitor_reactions))    #just set to 0.3 since we do not have MC simulations


    # Set up lists to hold output data
    output_foil_index = []
    output_mu = []
    output_unc_mu = []
    output_percent_unc = []


    # Get correlation submtarix for each monitor foil - n x n, where n= # of reactions per foil
    for i in range(0, number_of_monitor_foils): # Monitor reactions per foil [3,3,1]
        submatrix = np.ones((monitor_reactions_per_foil[i], monitor_reactions_per_foil[i]))
        corr_areal_density[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = submatrix
        corr_reaction_integral[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = 0.3*submatrix


    # Ensure all diagonal elements are still ones
    np.fill_diagonal(corr_lambda,1)
    np.fill_diagonal(corr_areal_density,1)
    np.fill_diagonal(corr_reaction_integral,1)
    np.fill_diagonal(corr_t_irradiation,1)
    np.fill_diagonal(corr_EoB_activities,1)
    # print("corr_lambda")
    # print(corr_lambda)
    # print("corr_areal_density")
    # print(corr_areal_density)
    # print("corr_reaction_integral")
    # print(corr_reaction_integral)
    # print("corr_t_irradiation")
    # print(corr_t_irradiation)
    # print("corr_EoB_activities")
    # print(corr_EoB_activities)


    # Loop over all beam positions
    number_of_energies = len(areal_density)
    # print(number_of_energies)
    # Test mode!!!!
    # number_of_energies = 1

    # Hold curents as we go along...
    currents = np.zeros((number_of_energies,number_of_monitor_reactions))
    unc_currents = np.zeros((number_of_energies,number_of_monitor_reactions))
    # function_dictionary = {'dIdA0':dIdA0, 'dIdRhoDr':dIdRhoDr, 'dIdLambda':dIdLambda, 'dIdTIrradiation':dIdTIrradiation, 'dIdIntegral':dIdIntegral}
    function_dictionary = {'0':dIdA0, '1':dIdRhoDr, '2':dIdLambda, '3':dIdTIrradiation, '4':dIdIntegral}


    # print('ad: ',areal_density)
    # print('unc_ad: ',uncertainty_areal_density)
    # print('A0:',EoB_activities)
    # print('unc_A0: ',uncertainty_EoB_activities)
    # print('rxn_int: ',reaction_integral)
    # print('unc_rxn_int: ',unc_rxn_int)
    # print('delta_t: ',t_irradiation)
    # print('unc_delta_t: ',uncertainty_t_irradiation)
    # print('loop_lambdas: ',lambdas)
    # print('uncertainty_lambdas: ',uncertainty_lambdas)

    for i_energy in range(0, number_of_energies):
        #print('i_energy: ',i_energy)
        # Get nonzero entries in A0:
        nonzero_indices = np.nonzero(EoB_activities[i_energy,:])
        ad = areal_density[i_energy,:]
        unc_ad = uncertainty_areal_density[i_energy,:]
        A0 = EoB_activities[i_energy,:]
        unc_A0 = uncertainty_EoB_activities[i_energy,:]
        rxn_int = reaction_integral[i_energy,:]
        delta_t = np.ones(number_of_monitor_reactions) *t_irradiation[i_energy]
        # delta_t = t_irradiation[i_energy]
        unc_delta_t = np.ones(number_of_monitor_reactions) *uncertainty_t_irradiation[i_energy]
        # unc_delta_t = uncertainty_t_irradiation[i_energy]
        #percent_rn_uncertainties = np.array([0.051054188386, 0.064100768909, 0.084213661384, 0.04385708308])
        #unc_rxn_int = rxn_int *    percent_rn_uncertainties      #rxn = reactions
        unc_rxn_int = unc_rxn_integral[i_energy,:]
        loop_lambdas = lambdas
        uncertainty_lambdas = loop_lambdas * 0.001

        # print('ad: ',ad)
        # print('unc_ad: ',unc_ad)
        # print('A0:',A0)
        # print('unc_A0: ',unc_A0)
        # print('rxn_int: ',rxn_int)
        # print('unc_rxn_int: ',unc_rxn_int)
        # print('delta_t: ',delta_t)
        # print('unc_delta_t: ',unc_delta_t)
        # print('loop_lambdas: ',loop_lambdas)
        # print('uncertainty_lambdas: ',uncertainty_lambdas)

        if len(np.transpose(nonzero_indices)) == number_of_monitor_reactions:

            # No nonzero indices!!!
            # Keep normal correlation matrices
            loop_corr_lambda = corr_lambda
            loop_corr_areal_density = corr_areal_density
            loop_corr_reaction_integral = corr_reaction_integral
            loop_corr_t_irradiation = corr_t_irradiation
            loop_corr_EoB_activities = corr_EoB_activities

        else:
            # Some nonzero indices
            # Find which indices are missing!
            temp3 = np.asarray(nonzero_indices[0])
            temp4 = np.array(range(0, number_of_monitor_reactions))
            disjoint_indices = np.setdiff1d(temp4,temp3,assume_unique=False).tolist()
            # print('disjoint indices: ', disjoint_indices)

            # Delete rows and columns in correlation matries of disjoint indices
            # gen = (x for x in xyz if x not in a)
            if len(disjoint_indices) != 1:
                loop_corr_lambda = np.delete(corr_lambda,np.array(disjoint_indices),0)
                loop_corr_lambda = np.delete(loop_corr_lambda,np.array(disjoint_indices),1)
                loop_corr_areal_density = np.delete(corr_areal_density,np.array(disjoint_indices),0)
                loop_corr_areal_density = np.delete(loop_corr_areal_density,np.array(disjoint_indices),1)
                loop_corr_reaction_integral = np.delete(corr_reaction_integral,np.array(disjoint_indices),0)
                loop_corr_reaction_integral = np.delete(loop_corr_reaction_integral,np.array(disjoint_indices),1)
                loop_corr_t_irradiation = np.delete(corr_t_irradiation,np.array(disjoint_indices),0)
                loop_corr_t_irradiation = np.delete(loop_corr_t_irradiation,np.array(disjoint_indices),1)
                loop_corr_EoB_activities = np.delete(corr_EoB_activities,np.array(disjoint_indices),0)
                loop_corr_EoB_activities = np.delete(loop_corr_EoB_activities,np.array(disjoint_indices),1)
            else:
                for disjoint_index in disjoint_indices:
                    loop_corr_lambda = np.delete(corr_lambda,disjoint_index,0)
                    loop_corr_lambda = np.delete(loop_corr_lambda,disjoint_index,1)
                    loop_corr_areal_density = np.delete(corr_areal_density,disjoint_index,0)
                    loop_corr_areal_density = np.delete(loop_corr_areal_density,disjoint_index,1)
                    loop_corr_reaction_integral = np.delete(corr_reaction_integral,disjoint_index,0)
                    loop_corr_reaction_integral = np.delete(loop_corr_reaction_integral,disjoint_index,1)
                    loop_corr_t_irradiation = np.delete(corr_t_irradiation,disjoint_index,0)
                    loop_corr_t_irradiation = np.delete(loop_corr_t_irradiation,disjoint_index,1)
                    loop_corr_EoB_activities = np.delete(corr_EoB_activities,disjoint_index,0)
                    loop_corr_EoB_activities = np.delete(loop_corr_EoB_activities,disjoint_index,1)

        # print('beam_current inputs: ', A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int)
        temp_currents =  beam_current(A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int)
        currents[i_energy, :] =  temp_currents
        # print('temp_currents: ', temp_currents)
        # print('unc_beam_current inputs: ', A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int, unc_A0, unc_ad, uncertainty_lambdas[i_energy,:], unc_delta_t, unc_rxn_int)
        unc_temp_currents = sigma_I_approximate(A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int, unc_A0, unc_ad, uncertainty_lambdas[i_energy,:], unc_delta_t, unc_rxn_int)
        unc_currents[i_energy,:] = unc_temp_currents

        value_array = np.array([A0, ad, loop_lambdas[0], delta_t, rxn_int])
        uncertainty_array = np.array([unc_A0, unc_ad, uncertainty_lambdas[0], unc_delta_t, unc_rxn_int])
        correlation_array = np.array([loop_corr_EoB_activities, loop_corr_areal_density, loop_corr_lambda, loop_corr_t_irradiation, loop_corr_reaction_integral])


        # Set up covariance matrix for current energy position
        cov = np.zeros((len(nonzero_indices[0]),len(nonzero_indices[0])))

        # NaN handling - replace range(0,number_of_monitor_reactions) with indices of nonzero elements of A0?
        # Fill correlation matrices
        for i_index,i_element in enumerate(nonzero_indices[0]):
        # for i in range(0,number_of_monitor_reactions):
            for j_index,j_element in enumerate(nonzero_indices[0]):
            # for j in range(0, number_of_monitor_reactions):
                for dict_index,dict_key in enumerate(function_dictionary):
                    # print("dict_index: ",dict_index)
                    # print("dict_value: ",function_dictionary[dict_key])
                    # print(type(dict_key))
                    # print(A0[i_element], ad[i_element], loop_lambdas[0,i_element], delta_t[i_element], rxn_int[i_element])
                    dIdxi = function_dictionary[dict_key](A0[i_element], ad[i_element], loop_lambdas[0,i_element], delta_t[i_element], rxn_int[i_element])
                    dIdxj = function_dictionary[dict_key](A0[j_element], ad[j_element], loop_lambdas[0,j_element], delta_t[j_element], rxn_int[j_element])
                    # print("dIdx_i: ",dIdxi)
                    # print("dIdx_j: ",dIdxj)
                    # print("unc_xi: ",uncertainty_array[dict_index,i_element])
                    # print("unc_xj: ",uncertainty_array[dict_index,j_element])
                    # print("corr_x: ", correlation_array[dict_index,i_index,j_index])
                    cov[i_index,j_index] += dIdxi * uncertainty_array[int(dict_key),i_element] *  correlation_array[int(dict_key),i_index,j_index] *  uncertainty_array[int(dict_key),j_element] * dIdxj

        # print("Final covariance matrix: \n", cov)
        inverted_covariance = np.linalg.inv(cov)
        numerator = 0.0
        denominator = 0.0

        for i_index,i_element in enumerate(nonzero_indices[0]):
            for j_index,j_element in enumerate(nonzero_indices[0]):
                numerator += temp_currents[j_element] * inverted_covariance[i_index,j_index]
                denominator +=  inverted_covariance[i_index,j_index]

        weighted_average_current = numerator/denominator
        uncertainty_weighted_average_current = np.sqrt(1.0/denominator)


        ####print("weighted_average_current: ",weighted_average_current, " +/- ",uncertainty_weighted_average_current, " nA     (", 100*uncertainty_weighted_average_current/weighted_average_current ," %)")

        # Append values for current energy
        output_foil_index.append(i_energy)
        output_mu.append(weighted_average_current)
        output_unc_mu.append(uncertainty_weighted_average_current)
        output_percent_unc.append(100*uncertainty_weighted_average_current/weighted_average_current)
        ####print("********************************************************************\n")
    # print("Raw currents: \n",currents)
    # print("Raw unc_currents: \n",unc_currents)


    #matlab_avg_currents = np.array(  [103.6055,   99.3354,  104.7844,  109.5334,  103.4154,   97.9257,   92.6829,   90.6173,   92.6035,   88.0497,   87.8754,   75.7701,   66.8488,    0.1905])
    #matlab_unc_avg_currents = np.array([3.3891,    2.9333,    3.3763,    3.4887,     3.3309,     3.2648,     3.8734,     2.8113,     2.8976,     3.3152,     3.6999,     2.3559,     2.6791,     0.2861])

    # Save final results to csv
    outfile = np.stack((np.transpose(output_foil_index),np.transpose(output_mu),np.transpose(output_unc_mu),np.transpose(output_percent_unc)), axis=-1)
    

    

    #import os
    #path = os.getcwd()
    #print("Does it save? ")
    #print(csv_filename)
    #save_string = "weighted_BC_"+csv_filename
    #print(save_string)

    csv_outname = 'WABC_' + csv_filename[10:-11] + '.csv'
    
    if save_csv==True:
        #np.savetxt("./{}".format(csv_filename), outfile, delimiter=",", header="Foil Index, Average Current (nA), Uncertainty in Average Current (nA), % Uncertainty")
        np.savetxt("./{}".format(csv_outname), outfile, delimiter=",", header="Foil Index, Average Current (nA), Uncertainty in Average Current (nA), % Uncertainty")
        #np.savetxt("weighted_BC_{}".format(csv_filename), outfile, delimiter=",", header="Foil Index, Average Current (nA), Uncertainty in Average Current (nA), % Uncertainty")



    # Plot output comparisons
    #
    # plt.gca().set_prop_cycle(None)
    output_foil_index2 = np.array(output_foil_index) - 0.2

    # plt.clf()

    # plt.gca()
    plt.errorbar(output_foil_index, output_mu, yerr=output_unc_mu, capsize=10.0, markersize=4.0,  marker='s', ls=' ', color='black', capthick=1.5, linewidth=3.0)
    # plt.errorbar(output_foil_index2, matlab_avg_currents, yerr=matlab_unc_avg_currents, capsize=10.0, capthick=2.0, markersize=8.0,  marker='.', ls=' ',  linewidth=2.0)
    for i in range(len(output_foil_index)):
        plt.errorbar(np.ones(len(currents[i,:]))*(0.2+output_foil_index[i]), currents[i,:], markersize=4.0,  yerr=unc_currents[i,:], capsize=5.0,  marker='.', ls=' ',linewidth=0.5, capthick=0.5, color='red')
    plt.legend(['Average Currents', 'Individual Currents'],loc='lower left')
    # plt.legend(['Actual Currents', 'Approximate Currents', 'Individual Currents'],loc='lower left')
    plt.show()
    

    #output_mu = output_mu.reverse()

    return output_mu[::-1], output_unc_mu[::-1 ] #returning reversed lists







curr, unc_curr = Average_BeamCurrent(A0, sigma_A0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time, csv_filename='averaged_currents.csv', save_csv=False)


print(curr)
print(unc_curr)
