import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge
from statistics import stdev

#________________________________________________________________________30 MeV ______________________________________________________________________________________________
number_of_monitor_foils_30MeV = 2
monitor_reactions_per_foil_30MeV = np.array([3,2]) 
number_of_comp_30MeV = 5
number_of_monitor_reactions_30MeV = 5 

A0_30MeV = np.array([ # [Bq]
    [261.812792566163,   3735.718330951421,  205620.8434783987, 301.4015855517144,  12360.35920909536],    # Comp 1 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
    [195.41965991684066, 5708.687427808714,  196653.5010576076, 316.4517638978664,   21399.665055428068],    # Comp 2 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
    [507.91693209926547, 4671.007978127567,  326995.26336984494, 423.3649475924, 10083.278760252526],    # Comp 3 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
    [642.8379460389074,  1094.8517316309733, 654013.0791957023,  20.57305322754311, 600.0895099486978]     # Comp 4 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
sigma_A0_30MeV = np.array([ # [Bq]
    [17.582548417585322, 206.4249733635162,  4144.3555481768835,  3.4339708324821996,  81.66495061632561],      # Comp 1 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1.9411697079227395, 831.5006106743914,  3649.5885250219026,  3.911924584714338, 259.31980441421035],    # Comp 2 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [8.772734507668872,  9.266332716923932,  5303.582226233545,  7.011882019032889,  156.63991477782636],    # Comp 3 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [15.030405986416387, 402.34839028951893, 18690.088639284084, 0.6210899158274328, 6.903173653751662]      # Comp 4 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
])  

mass_density_30MeV = np.array([ # [nuclei/cm^2]                    
    [2.385836227282795e+20,  2.385836227282795e+20,  2.385836227282795e+20,  1.4512167908580024e+20, 1.4512167908580024e+20], # Comp 1 mass dens for Ni01, Ni01, Ni01, Ti01, Ti01
    [2.3704457213656056e+20, 2.3704457213656056e+20, 2.3704457213656056e+20, 1.4770078131489335e+20, 1.4770078131489335e+20], # Comp 2 mass dens for Ni02, Ni02, Ni02, Ti02, Ti02
    [2.366444189827136e+20,  2.366444189827136e+20,  2.366444189827136e+20,  1.418254606271544e+20,  1.418254606271544e+20],  # Comp 3 mass dens for Ni03, Ni03, Ni03, Ti03, Ti03
    [2.380500851898169e+20,  2.380500851898169e+20,  2.380500851898169e+20,  1.396363689790461e+20,  1.396363689790461e+20]  # Comp 4 mass dens for Ni04, Ni04, Ni04, Ti04, Ti04
])  
  
sigma_mass_density_30MeV = np.array([ # [nuclei/cm^2]  
    [1.5389502640420342e+17, 1.5389502640420342e+17, 1.5389502640420342e+17, 1.1952259940967288e+18, 1.1952259940967288e+18], # Comp 1 mass dens unc for Ni01, Ni01, Ni01, Ti01, Ti01
    [4.0013449889007386e+17, 4.0013449889007386e+17, 4.0013449889007386e+17, 1.1448329143295144e+18, 1.1448329143295144e+18], # Comp 2 mass dens unc for Ni02, Ni02, Ni02, Ti02, Ti02
    [7.080419383273164e+17,  7.080419383273164e+17,  7.080419383273164e+17,  1.081990495922333e+18,  1.081990495922333e+18],  # Comp 3 mass dens unc for Ni03, Ni03, Ni03, Ti03, Ti03
    [7.079628121784579e+17,  7.079628121784579e+17,  7.079628121784579e+17,  1.1574495386179843e+18, 1.1574495386179843e+18] # Comp 4 mass dens unc for Ni04, Ni04, Ni04, Ti04, Ti04
])   
  
lambda__30MeV = np.log(2)/np.array([ # [1/s]                                        # each row will be the same
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 1 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 2 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 3 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600]     # Comp 4 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
reaction_integral_30MeV = 1e4*np.array([ # [cm^2]
    [9.887459938071693e-31, 1.5345064286646657e-29, 1.4221729159539438e-30,  2.419567936368279e-30, 1.8228629152784894e-29],  # Comp 1 int for 56Co, 58Co, 61Cu, 46Sc, 48V
    [7.512261908970988e-31, 2.1537694649567716e-29,  1.7287059859684096e-30,  2.555534220184229e-30, 3.1984734875795636e-29],  # Comp 2 int for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1.9696922400525485e-30, 1.7946307709654718e-29, 2.58334360353661e-30, 3.317125913498328e-30,  1.565598289435036e-29],  # Comp 3 int for 56Co, 58Co, 61Cu, 46Sc, 48V
    [3.058762296497243e-30, 4.636953104138198e-30,  6.2163585870201584e-30,  3.353434333694576e-31, 9.071677145107408e-31]   # Comp 4 int for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
uncertainty_integral_30MeV = 1e4*np.array([ # [cm^2]
    [5.58566107107498e-32,  8.286376381070242e-31,  7.597706274519804e-32,  9.985248787957728e-32,  9.43958989512594e-31],  # Comp 1 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [4.67681957513234e-32, 1.1614983555781079e-30,  8.947601283325482e-32,  1.052362574869065e-31, 1.6516910215469435e-30],  # Comp 2 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1.159623146291675e-31, 1.0755624509537944e-30,  1.2308573518244383e-31, 1.4078351991377597e-31, 8.383723215823517e-31],  # Comp 3 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1.826727140117606e-31, 3.711658408563801e-31,  3.541362079215962e-31, 1.938911406538652e-32, 9.702245157980755e-32]  # Comp 4 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
irr_time_30MeV = np.ones((number_of_comp_30MeV, number_of_monitor_reactions_30MeV))*1200 # [s]
sigma_irr_time_30MeV = np.ones((number_of_comp_30MeV, number_of_monitor_reactions_30MeV))*3 # [s]








#________________________________________________________________________50 MeV ______________________________________________________________________________________________

number_of_monitor_foils_50MeV = 2
monitor_reactions_per_foil_50MeV = np.array([1,2])
number_of_comp_50MeV = 5
number_of_monitor_reactions_50MeV = 3


A0_50MeV = np.array([ # [Bq]
    [660.6582922597377,  682.303263776218,   2676.4044753380927],     # Comp 1 activities for 56Co, 46Sc, 48V
    [867.2472346205243,  446.4477144386634,  3515.5863153648215],    # Comp 2 activities for 56Co, 46Sc, 48V
    [1101.694341853326,  305.2401745402718,  4435.525810342462],     # Comp 3 activities for 56Co, 46Sc, 48V
    [1404.4599460617442, 238.7734590693519,  5943.619406742724],     # Comp 4 activities for 56Co, 46Sc, 48V
    [2169.269618542412,  214.6928704035015,  9912.983408250037]      # Comp 5 activities for 56Co, 46Sc, 48V
])  
  
sigma_A0_50MeV = np.array([ # [Bq]
    [46.66272559761359,  10.405047778884516,  27.816502688756184],   # Comp 1 activities unc for 56Co, 46Sc, 48V
    [62.05286622840644,  2.8470515224086257,  27.614757458099724],    # Comp 2 activities unc for 56Co, 46Sc, 48V
    [78.89663235383438,  0.3523340356073281,  30.840055999870838],   # Comp 3 activities unc for 56Co, 46Sc, 48V
    [98.24060397528628,  1.672518728554704,   54.22565920058882],    # Comp 4 activities unc for 56Co, 46Sc, 48V
    [150.1916177435665,  2.229527535870038,   19.108166277136167]    # Comp 5 activities unc for 56Co, 46Sc, 48V
])  

mass_density_50MeV = np.array([ # [nuclei/cm^2]                      
    [2.163312142336825e+20,  1.3826504145235758e+20, 1.3826504145235758e+20],    # Comp 1 mass density for Fe01, Ti06, Ti06
    [2.174095812853434e+20,  1.4040380915453236e+20, 1.4040380915453236e+20],    # Comp 2 mass density for Fe02, Ti08, Ti08
    [2.1768995671877517e+20, 1.4133480215430254e+20, 1.4133480215430254e+20],    # Comp 3 mass density for Fe03, Ti09, Ti09
    [2.1675177738383022e+20, 1.4139770708671947e+20, 1.4139770708671947e+20],    # Comp 4 mass density for Fe04, Ti10, Ti10
    [2.1683804674796313e+20, 1.3840343230367475e+20, 1.3840343230367475e+20]     # Comp 5 mass density for Fe05, Ti11, Ti11
])  
  
sigma_mass_density_50MeV = np.array([ # [nuclei/cm^2]  
    [3.2350795013732275e+17, 1.0693457318535995e+18, 1.0693457318535995e+18],    # Comp 1 mass density unc for Fe01, Ti06, Ti06
    [3.344665770589681e+17,  1.308426385320377e+18,  1.308426385320377e+18],     # Comp 2 mass density unc for Fe02, Ti08, Ti08
    [1.2077690426674125e+18, 1.0819219395230024e+18, 1.0819219395230024e+18],    # Comp 3 mass density unc for Fe03, Ti09, Ti09
    [5.3933428859724294e+17, 3.0319924720390205e+18, 3.0319924720390205e+18],    # Comp 4 mass density unc for Fe04, Ti10, Ti10
    [1.2078128857102303e+18, 1.8880442303825955e+17, 1.8880442303825955e+17]     # Comp 5 mass density unc for Fe05, Ti11, Ti11
])  

  
lambda__50MeV = np.log(2)/np.array([ # [1/s]    # each row will be the same
    [77.236*24*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 1 lambda for 56Co, 46Sc, 48V
    [77.236*24*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 2 lambda for 56Co, 46Sc, 48V
    [77.236*24*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 3 lambda for 56Co, 46Sc, 48V
    [77.236*24*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 4 lambda for 56Co, 46Sc, 48V
    [77.236*24*3600, 83.79*24*3600, 15.9735*24*3600]     # Comp 5 lambda for 56Co, 46Sc, 48V
])  
  
reaction_integral_50MeV = 1e4*np.array([ # [cm^2]
    [4.2801999070724114e-30, 6.940858820023156e-30,  5.609459687348817e-30], # Comp 1 int for 56Co, 46Sc, 48V
    [5.353658424185412e-30,  4.7202060602526754e-30,  7.2559378069068e-30],  # Comp 2 int for 56Co, 46Sc, 48V
    [6.742712662562115e-30,  3.1113272630680595e-30, 8.808775181686762e-30],  # Comp 3 int for 56Co, 46Sc, 48V
    [9.589102219731265e-30,  2.5162052569096332e-30, 1.215999082863533e-29], # Comp 4 int for 56Co, 46Sc, 48V
    [1.659191705527795e-29,  2.386390265266877e-30,  2.1262639306973067e-29]   # Comp 5 int for 56Co, 46Sc, 48V
])  
  
uncertainty_integral_50MeV = 1e4*np.array([ # [cm^2]
    [4.172834386347752e-31,  3.1458335207044577e-31, 3.8271924485962245e-31],  # Comp 1 int unc for 56Co, 46Sc, 48V
    [4.912465966010479e-31,  2.0451997740209564e-31,  4.003426777795037e-31],  # Comp 2 int unc for 56Co, 46Sc, 48V
    [5.797528842937585e-31,  1.3178490858385847e-31, 4.683974407660781e-31],   # Comp 3 int unc for 56Co, 46Sc, 48V
    [7.174470959280437e-31,  1.0476196619148873e-31,  6.387021968198101e-31],   # Comp 4 int unc for 56Co, 46Sc, 48V
    [1.0174139999477894e-30, 9.835068351272783e-32,  1.100820065367547e-30]   # Comp 5 int unc for 56Co, 46Sc, 48V
])   
  
irr_time_50MeV = np.ones((number_of_comp_50MeV, number_of_monitor_reactions_50MeV))*1200 # [s]
sigma_irr_time_50MeV = np.ones((number_of_comp_50MeV, number_of_monitor_reactions_50MeV))*3 # [s]

#_____________________________________________________________________________________________________________________________________________________________________________






# A0 = A0_30MeV
# sigma_A0 = sigma_A0_30MeV 
# mass_density = mass_density_30MeV
# sigma_mass_density = sigma_mass_density_30MeV 
# lambda_= lambda__30MeV 
# reaction_integral = reaction_integral_30MeV
# uncertainty_integral = uncertainty_integral_30MeV 
# irr_time = irr_time_30MeV 
# sigma_irr_time = sigma_irr_time_30MeV
# number_of_monitor_foils = number_of_monitor_foils_30MeV
# monitor_reactions_per_foil = monitor_reactions_per_foil_30MeV


A0 = A0_50MeV
sigma_A0 = sigma_A0_50MeV 
mass_density = mass_density_50MeV
sigma_mass_density = sigma_mass_density_50MeV 
lambda_= lambda__50MeV 
reaction_integral = reaction_integral_50MeV
uncertainty_integral = uncertainty_integral_50MeV 
irr_time = irr_time_50MeV 
sigma_irr_time = sigma_irr_time_50MeV
number_of_monitor_foils = number_of_monitor_foils_50MeV
monitor_reactions_per_foil = monitor_reactions_per_foil_50MeV










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
    #unc_rxn_integral = uncertainty_integral #rxn_int * 	percent_rn_uncertainties      #rxn = reactions
    unc_rxn_integral = uncertainty_integral * reaction_integral#rxn_int * 	percent_rn_uncertainties      #rxn = reactions
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
		#unc_rxn_int = rxn_int * 	percent_rn_uncertainties      #rxn = reactions
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
    # plt.ylim(-5,170)
	# plt.legend(['Actual Currents', 'Approximate Currents', 'Individual Currents'],loc='lower left')
    plt.show()
    

    #output_mu = output_mu.reverse()

    return output_mu, output_unc_mu







curr, unc_curr = Average_BeamCurrent(A0, sigma_A0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time, csv_filename='averaged_currents.csv', save_csv=False)

print('before adding std to unc:')
print('bc: ', curr)
print('unc_bc: ', unc_curr)





# # Adding sqrt(std) to the unc for 30 MeV stack:


# bc_dict = {'56CO': {'beam_current': [142.6690228199885, 141.07493375244087, 140.07551879558267, 113.488439539866], 
#                     'beam_current_unc': [16.13653401452156, 17.586330818269932, 16.52981786203812, 13.584181291586818], 
#                     'energy': [27.35133924, 20.450330565, 14.083908879999997, 8.848522071305181]}, 
#            '58CO': {'beam_current': [120.34076867910537, 131.872740763852, 129.71411332064656, 116.97748057168313, 468.31173326214855], 
#                     'beam_current_unc': [13.01330924777981, 14.247413889265868, 15.582529477671539, 18.750063129542692, 108.46808283879223], 
#                     'energy': [27.35133924, 20.450330565, 14.083908879999997, 8.848522071305181, 2.3286219274448703]}, 
#            '61CU': {'beam_current': [145.22318046876342, 115.00384784884434, 128.18135577071817, 105.91168407157757, 6.009329520490208], 
#                     'beam_current_unc': [15.54808323491192, 11.936098772249416, 12.268685839683723, 12.10442484188049, 1.3459110723307477], 
#                     'energy': [27.35133924, 20.450330565, 14.083908879999997, 8.848522071305181, 2.3286219274448703]}, 
#            '48V': {'beam_current': [124.24905346888342, 120.45635827792036, 120.75789499180249, 125.97321213089911, 102.47112196443523], 
#                     'beam_current_unc': [13.044917518433676, 12.594543763041111, 13.077666870850727, 27.03408006522823, 25.45551112893367], 
#                     'energy': [25.542917359999997, 18.134367400000002, 10.866410484000001, 3.9274296448168275, 1.3297790055248617]}, 
#            '46SC': {'beam_current': [119.7043426397923, 116.91696038577017, 125.49696608010343, 61.26957739268947, 290.8286891311535], 
#                     'beam_current_unc': [10.09333188333647, 9.816366726559183, 10.842103608803457, 7.164281333735599, 139.05734096258246], 
#                     'energy': [25.542917359999997, 18.134367400000002, 10.866410484000001, 3.9274296448168275, 1.3297790055248617]}}



# comp1_list = []
# comp2_list = []
# comp3_list = []
# comp4_list = []


# for reaction in bc_dict:
#     comp1_list.append(bc_dict[reaction]['beam_current'][0])
#     comp2_list.append(bc_dict[reaction]['beam_current'][1])
#     comp3_list.append(bc_dict[reaction]['beam_current'][2])
#     comp4_list.append(bc_dict[reaction]['beam_current'][3])


# std1 = np.std(comp1_list)
# std2 = np.std(comp2_list)
# std3 = np.std(comp3_list)
# std4 = np.std(comp4_list)
# std_list = []

# std_list.append(std1)
# std_list.append(std2)
# std_list.append(std3)
# std_list.append(std4)

# final_unc_curr = [np.sqrt(std_dev + other_unc**2) for std_dev, other_unc in zip(std_list, unc_curr)]


# print('after adding std to unc:')
# print('bc: ', curr)
# print('unc_bc: ', final_unc_curr)




#Adding sqrt(std) to the unc for 50 MeV stack:

bc_dict = {'56CO': {'beam_current': [91.71868901291089, 95.78064880345202, 96.48332515195193, 86.86262139072429, 77.50784615270966], 
                    'beam_current_unc': [17.891687001061534, 17.586591409360715, 16.633314233940332, 13.01251633646137, 9.55268221176066], 
                    'energy': [48.2318781, 41.779047799999994, 37.08828046666665, 31.872687066666668, 25.879117466666663]}, 
           '46SC': {'beam_current': [99.14854343832273, 93.94317721163227, 96.80134065374573, 93.59057325769712, 90.64915956974917], 
                    'beam_current_unc': [9.131366995365276, 8.340736045381597, 8.347731088799447, 8.77906515572666, 7.490193553873115], 
                    'energy': [47.09795026666667, 40.50592943333334, 35.68120613333333, 30.28116373333334, 23.99016156666667]}, 
           '48V': {'beam_current': [91.76279410520966, 91.76427218744738, 94.73899437338731, 91.92289189263836, 89.57537766869422], 
                    'beam_current_unc': [12.610055941263163, 10.279792929471318, 10.190216978411376, 10.440332545328687, 9.289141582970178], 
                    'energy': [47.09795026666667, 40.50592943333334, 35.68120613333333, 30.28116373333334, 23.99016156666667]}}


comp1_list = []
comp2_list = []
comp3_list = []
comp4_list = []
comp5_list = []



for reaction in bc_dict:
    comp1_list.append(bc_dict[reaction]['beam_current'][0])
    comp2_list.append(bc_dict[reaction]['beam_current'][1])
    comp3_list.append(bc_dict[reaction]['beam_current'][2])
    comp4_list.append(bc_dict[reaction]['beam_current'][3])
    comp5_list.append(bc_dict[reaction]['beam_current'][4])



std1 = np.std(comp1_list)
std2 = np.std(comp2_list)
std3 = np.std(comp3_list)
std4 = np.std(comp4_list)
std5 = np.std(comp5_list)
std_list = []

std_list.append(std1)
std_list.append(std2)
std_list.append(std3)
std_list.append(std4)
std_list.append(std5)


final_unc_curr = [np.sqrt(std_dev + other_unc**2) for std_dev, other_unc in zip(std_list, unc_curr)]


print('after adding std to unc:')
print('bc: ', curr)
print('unc_bc: ', final_unc_curr)

