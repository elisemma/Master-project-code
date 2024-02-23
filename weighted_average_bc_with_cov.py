import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge

#________________________________________________________________________30 MeV ______________________________________________________________________________________________
number_of_monitor_foils_30MeV = 2
monitor_reactions_per_foil_30MeV = np.array([3,2]) 
number_of_comp_30MeV = 5
number_of_monitor_reactions_30MeV = 5 

A0_30MeV = np.array([ # [Bq]
    [261.812792566163,   3736.2259003636364, 220380.49108977473, 301.4362170806085,  12367.810902410458],    # Comp 1 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
    [195.41965991684066, 5709.463064094319,  210769.46467086408, 316.488124679544,   21412.566274580688],    # Comp 2 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
    [507.91693209926547, 4671.642624718161,  350467.27527558757, 423.41359286337206, 10089.357669794123],    # Comp 3 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
    [642.8379460389074,  1095.00048859536,   700958.7214614445,  20.575417102354017, 600.4512860881874],     # Comp 4 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1e-50,                  25.452633286057697*12.642, 461.7853246484494*12.642,  1e-50, 1e-50]     # Comp 5 activities for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
sigma_A0_30MeV = np.array([ # [Bq]
    [17.582548417585322, 206.45302017891137, 4441.841081947311,  3.434365403830543,  81.7141847549622],      # Comp 1 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1.9411697079227395, 831.6135862144428,  3911.559231355066,  3.9123740720603895, 259.47614231008066],    # Comp 2 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [8.772734507668872,  9.267591691987594,  5684.278089150438,  7.012687689945556,  156.73434909843314],    # Comp 3 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [15.030405986416387, 402.40305716671895, 20031.679842941638, 0.6211612800249608, 6.907335400268808],     # Comp 4 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1e12,                  2.465479823966241*12.642,  118.4114034137288*12.642,  1e12, 1e12]     # Comp 5 activities unc for 56Co, 58Co, 61Cu, 46Sc, 48V
])  

# mass_density_30MeV = np.array([                       
#     [23.253, 23.253, 23.253, 11.535, 11.535],    # Comp 1 mass density for Ni01, Ni01, Ni01, Ti01, Ti01
#     [23.103, 23.103, 23.103, 11.740, 11.740],    # Comp 2 mass density for Ni02, Ni02, Ni02, Ti02, Ti02
#     [23.064, 23.064, 23.064, 11.273, 11.273],    # Comp 3 mass density for Ni03, Ni03, Ni03, Ti03, Ti03
#     [23.201, 23.201, 23.201, 11.099, 11.099],    # Comp 4 mass density for Ni04, Ni04, Ni04, Ti04, Ti04
#     [22.746, 22.746, 22.746, 11.317, 11.317]     # Comp 5 mass density for Ni05, Ni05, Ni05, Ti05, Ti05
# ])  
  
# sigma_mass_density_30MeV = np.array([
#     [0.015, 0.015, 0.015, 0.095, 0.095],         # Comp 1 mass density unc for Ni01, Ni01, Ni01, Ti01, Ti01
#     [0.039, 0.039, 0.039, 0.091, 0.091],         # Comp 2 mass density unc for Ni02, Ni02, Ni02, Ti02, Ti02
#     [0.069, 0.069, 0.069, 0.086, 0.086],         # Comp 3 mass density unc for Ni03, Ni03, Ni03, Ti03, Ti03
#     [0.069, 0.069, 0.069, 0.092, 0.092],         # Comp 4 mass density unc for Ni04, Ni04, Ni04, Ti04, Ti04
#     [0.019, 0.019, 0.019, 0.030, 0.030]          # Comp 5 mass density unc for Ni05, Ni05, Ni05, Ti05, Ti05
# ]) 

mass_density_30MeV = np.array([ # [nuclei/cm^2]                    
    [2.385836227282795e+20,  2.385836227282795e+20,  2.385836227282795e+20,  1.4512167908580024e+20, 1.4512167908580024e+20], # Comp 1 mass dens for Ni01, Ni01, Ni01, Ti01, Ti01
    [2.3704457213656056e+20, 2.3704457213656056e+20, 2.3704457213656056e+20, 1.4770078131489335e+20, 1.4770078131489335e+20], # Comp 2 mass dens for Ni02, Ni02, Ni02, Ti02, Ti02
    [2.366444189827136e+20,  2.366444189827136e+20,  2.366444189827136e+20,  1.418254606271544e+20,  1.418254606271544e+20],  # Comp 3 mass dens for Ni03, Ni03, Ni03, Ti03, Ti03
    [2.380500851898169e+20,  2.380500851898169e+20,  2.380500851898169e+20,  1.396363689790461e+20,  1.396363689790461e+20],  # Comp 4 mass dens for Ni04, Ni04, Ni04, Ti04, Ti04
    [2.333816317282693e+20,  2.333816317282693e+20,  2.333816317282693e+20,  1.423790240324232e+20,  1.423790240324232e+20]   # Comp 5 mass dens for Ni05, Ni05, Ni05, Ti05, Ti05
])  
  
sigma_mass_density_30MeV = np.array([ # [nuclei/cm^2]  
    [1.5389502640420342e+17, 1.5389502640420342e+17, 1.5389502640420342e+17, 1.1952259940967288e+18, 1.1952259940967288e+18], # Comp 1 mass dens unc for Ni01, Ni01, Ni01, Ti01, Ti01
    [4.0013449889007386e+17, 4.0013449889007386e+17, 4.0013449889007386e+17, 1.1448329143295144e+18, 1.1448329143295144e+18], # Comp 2 mass dens unc for Ni02, Ni02, Ni02, Ti02, Ti02
    [7.080419383273164e+17,  7.080419383273164e+17,  7.080419383273164e+17,  1.081990495922333e+18,  1.081990495922333e+18],  # Comp 3 mass dens unc for Ni03, Ni03, Ni03, Ti03, Ti03
    [7.079628121784579e+17,  7.079628121784579e+17,  7.079628121784579e+17,  1.1574495386179843e+18, 1.1574495386179843e+18], # Comp 4 mass dens unc for Ni04, Ni04, Ni04, Ti04, Ti04
    [1.9488015308666957e+17, 1.9488015308666957e+17, 1.9488015308666957e+17, 3.774585126968663e+17,  3.774585126968663e+17]   # Comp 5 mass dens unc for Ni05, Ni05, Ni05, Ti05, Ti05
])   
  
lambda__30MeV = np.log(2)/np.array([ # [1/s]                                        # each row will be the same
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 1 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 2 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 3 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600],    # Comp 4 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
    [77.236*24*3600, 70.86*24*3600, 3.339*3600, 83.79*24*3600, 15.9735*24*3600]     # Comp 5 lambda for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
reaction_integral_30MeV = 1e4*np.array([ # [cm^2]
    [9.868468789066539e-31, 1.5361098661267076e-29, 1.422638698833147e-30,  2.4190166123842706e-30, 1.8271896436850462e-29],  # Comp 1 int for 56Co, 58Co, 61Cu, 46Sc, 48V
    [7.546457582261141e-31, 2.157347386230939e-29,  1.732476913736326e-30,  2.5656634319791498e-30, 3.2051513865015736e-29],  # Comp 2 int for 56Co, 58Co, 61Cu, 46Sc, 48V
    [2.0114277147898e-30,   1.7671033103200475e-29, 2.6169459389045837e-30, 3.273026664765459e-30,  1.4826991003748955e-29],  # Comp 3 int for 56Co, 58Co, 61Cu, 46Sc, 48V
    [3.006599245089043e-30, 4.299345718197932e-30,  6.351913050401064e-30,  2.8332737095796377e-31, 8.057275110831515e-31],   # Comp 4 int for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1,                     3.3735628091502424e-31, 8.16618313949607e-31,   1.203576743075432e-33,  1.9696939409830086e-33]   # Comp 5 int for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
uncertainty_integral_30MeV = 1e4*np.array([ # [cm^2]
    [5.574697852499285e-32,  8.295029704286848e-31,  7.600191560063854e-32,  9.981951030650341e-32,  9.462243749256337e-31],  # Comp 1 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [4.6950384114971025e-32, 1.164031694075335e-30,  8.959731837907958e-32,  1.0568544444325924e-31, 1.655751320482127e-30],  # Comp 2 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1.1823070709015115e-31, 1.0616938940527619e-30, 1.2453802324674207e-31, 1.3928406683925334e-31, 7.981419700271481e-31],  # Comp 3 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1.8007222096321475e-31, 3.480855851144396e-31,  3.7336235341067837e-31, 1.7276027933875856e-32, 8.584850901614173e-32],  # Comp 4 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
    [1,                      3.916078608773901e-32,  9.173104893466661e-32,  3.0849095349218334e-34, 2.8605518534780288e-34]  # Comp 5 int unc for 56Co, 58Co, 61Cu, 46Sc, 48V
])  
  
irr_time_30MeV = np.ones((number_of_comp_30MeV, number_of_monitor_reactions_30MeV))*1200 # [s]
sigma_irr_time_30MeV = np.ones((number_of_comp_30MeV, number_of_monitor_reactions_30MeV))*3 # [s]








#________________________________________________________________________50 MeV ______________________________________________________________________________________________

number_of_monitor_foils_50MeV = 2
monitor_reactions_per_foil_50MeV = np.array([1,2])
number_of_comp_50MeV = 5
number_of_monitor_reactions_50MeV = 3


A0_50MeV = np.array([ # [Bq]
    [655.9564121176588,  677.82608908275,    2585.555092750958],     # Comp 1 activities for 56Co, 46Sc, 48V
    [861.0750687694533,  443.5181894942165,  3396.2512716211227],    # Comp 2 activities for 56Co, 46Sc, 48V
    [1093.8536247846293, 303.23723292703966, 4284.963821779039],     # Comp 3 activities for 56Co, 46Sc, 48V
    [1394.464457609836,  237.20666237218154, 5741.865838977189],     # Comp 4 activities for 56Co, 46Sc, 48V
    [2153.831008503944,  213.28408701930022, 9576.491511308159]      # Comp 5 activities for 56Co, 46Sc, 48V
])  
  
sigma_A0_50MeV = np.array([ # [Bq]
    [46.33062873992905,  10.336771342401113,  26.872283424993523],   # Comp 1 activities unc for 56Co, 46Sc, 48V
    [61.611238320462014, 2.828369582318218,   26.67738647271755],    # Comp 2 activities unc for 56Co, 46Sc, 48V
    [78.33512799780877,  0.35002207130816204, 29.793203526154276],   # Comp 3 activities unc for 56Co, 46Sc, 48V
    [97.54142930298822,  1.6615439025852172,  52.38499295320197],    # Comp 4 activities unc for 56Co, 46Sc, 48V
    [149.1227096661337,  2.2148977036748057,  18.459547881969502]    # Comp 5 activities unc for 56Co, 46Sc, 48V
])  

# mass_density_50MeV = np.array([                       
#     [20.061, 10.990, 10.990],    # Comp 1 mass density for Fe01, Ti06, Ti06
#     [20.161, 11.160, 11.160],    # Comp 2 mass density for Fe02, Ti08, Ti08
#     [20.187, 11.234, 11.234],    # Comp 3 mass density for Fe03, Ti09, Ti09
#     [20.100, 11.239, 11.239],    # Comp 4 mass density for Fe04, Ti10, Ti10
#     [20.108, 11.001, 11.001]     # Comp 5 mass density for Fe05, Ti11, Ti11
# ])  
  
# sigma_mass_density_50MeV = np.array([
#     [0.036, 0.085, 0.085],         # Comp 1 mass density unc for Fe01, Ti06, Ti06
#     [0.031, 0.104, 0.104],         # Comp 2 mass density unc for Fe02, Ti08, Ti08
#     [0.112, 0.086, 0.086],         # Comp 3 mass density unc for Fe03, Ti09, Ti09
#     [0.050, 0.241, 0.241],         # Comp 4 mass density unc for Fe04, Ti10, Ti10
#     [0.112, 0.015, 0.015]          # Comp 5 mass density unc for Fe05, Ti11, Ti11
# ])  

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
    [4.2779983671409635e-30, 6.948483984428582e-30,  5.6063074104954894e-30], # Comp 1 int for 56Co, 46Sc, 48V
    [5.331767516943081e-30,  4.765629583471024e-30,  7.223346507166584e-30],  # Comp 2 int for 56Co, 46Sc, 48V
    [6.681272134443707e-30,  3.1537829296042605e-30, 8.734932501396683e-30],  # Comp 3 int for 56Co, 46Sc, 48V
    [9.408080530896652e-30,  2.5253018618753087e-30, 1.1918806006099854e-29], # Comp 4 int for 56Co, 46Sc, 48V
    [1.6007958861009106e-29, 2.394203550267864e-30,  2.047827429598293e-29]   # Comp 5 int for 56Co, 46Sc, 48V
])  
  
uncertainty_integral_50MeV = 1e4*np.array([ # [cm^2]
    [4.1685471808698985e-31, 3.1495005890023775e-31, 3.8293971120980756e-31],  # Comp 1 int unc for 56Co, 46Sc, 48V
    [4.897345721117305e-31,  2.066079066390724e-31,  3.9913225923279625e-31],  # Comp 2 int unc for 56Co, 46Sc, 48V
    [5.762174367947573e-31,  1.3369874371187185e-31, 4.646607238962349e-31],   # Comp 3 int unc for 56Co, 46Sc, 48V
    [7.095026475234454e-31,  1.052297613995038e-31,  6.271109628403489e-31],   # Comp 4 int unc for 56Co, 46Sc, 48V
    [9.932643683800369e-31,  9.869931597465798e-32,  1.0603634879519155e-30]   # Comp 5 int unc for 56Co, 46Sc, 48V
])   
  
irr_time_50MeV = np.ones((number_of_comp_50MeV, number_of_monitor_reactions_50MeV))*1200 # [s]
sigma_irr_time_50MeV = np.ones((number_of_comp_50MeV, number_of_monitor_reactions_50MeV))*3 # [s]

#_____________________________________________________________________________________________________________________________________________________________________________






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


# A0 = A0_50MeV
# sigma_A0 = sigma_A0_50MeV 
# mass_density = mass_density_50MeV
# sigma_mass_density = sigma_mass_density_50MeV 
# lambda_= lambda__50MeV 
# reaction_integral = reaction_integral_50MeV
# uncertainty_integral = uncertainty_integral_50MeV 
# irr_time = irr_time_50MeV 
# sigma_irr_time = sigma_irr_time_50MeV
# number_of_monitor_foils = number_of_monitor_foils_50MeV
# monitor_reactions_per_foil = monitor_reactions_per_foil_50MeV










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
    plt.ylim(-5,170)
	# plt.legend(['Actual Currents', 'Approximate Currents', 'Individual Currents'],loc='lower left')
    plt.show()
    

    #output_mu = output_mu.reverse()

    return output_mu[::-1], output_unc_mu[::-1 ] #returning reversed lists







curr, unc_curr = Average_BeamCurrent(A0, sigma_A0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time, csv_filename='averaged_currents.csv', save_csv=False)


print(curr)
print(unc_curr)
