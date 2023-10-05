import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


# csTi01 er bar en test. Automatisk kalibrering fungerer fint her
channel_array_csTi01 = np.array([1270, 2210, 2444, 2783, 5561])
energy_array_csTi01 = np.array([510.9989461, 889.277, 983.525, 1120.545, 2240.396])
energy_unc_array_csTi01 = np.array([0.0000031, 0.003, 0.004, 0.004, 0.01])
csTi01 =  "CS060317_Ti01_18cm_30MeV_calibration.Spe"


# Defining the channelnumber and real energy of the peaks I will use in the calibration
# for different foils
# channel_array_cjTi01 = np.array([1315, 2289, 2532, 2885, 5765])
# energy_array_cjTi01 = np.array([510.9989461, 889.277, 983.525, 1120.545, 2240.396])
# energy_unc_array_cjTi01 = np.array([0.0000031, 0.003, 0.004, 0.004, 0.01])
# cjTi01 = "CJ010317_Ti01_18cm_30MeV_calibration.Spe" 

channel_array_ckTi02 = np.array([1315, 2289, 2532, 2885, 3377, 5765]) #Using these values to make a calibration for all the Ti foils (CJ, CK, CL and CM)
energy_array_ckTi02 = np.array([510.9989461, 889.277, 983.525, 1120.545, 1312.105, 2240.396])
energy_unc_array_ckTi02 = np.array([0.0000031, 0.003, 0.004, 0.004, 0.004, 0.01])
subset_params_ckTi02 = np.array([50,50,20,50,70,50])
ckTi02 = "CK010317_Ti02_18cm_30MeV_calibration.Spe"

# channel_array_clTi03  = np.array([1315, 2289, 2885, 5765]) #Burde ha med en stor topp på ch = 3377
# energy_array_clTi03  = np.array([510.9989461, 889.277, 1120.545, 2240.396])
# energy_unc_array_clTi03  = np.array([0.0000031, 0.003, 0.004, 0.01])
# clTi03 = "CL010317_Ti03_18cm_30MeV_calibration.Spe"

# channel_array_cmTi04 = np.array([1315, 2289, 2532, 2885]) #Burde ha med en stor topp på ch = 3377
# energy_array_cmTi04 = np.array([510.9989461, 889.277, 983.525, 1120.545])
# energy_unc_array_cmTi04 = np.array([0.0000031, 0.003, 0.004, 0.004])
# cmTi04 = "CM010317_Ti04_18cm_30MeV_calibration.Spe"

channel_array_dfNi04 = np.array([1288, 2043, 2134, 6543]) #Using these values to makea calibration for all the Ni foils (DF and DG)
energy_array_dfNi04 = np.array([510.9989461, 810.7593, 846.770, 2598.500])
energy_unc_array_dfNi04 = np.array([0.0000031, 0.002, 0.002, 0.004])
subset_params_dfNi04 = np.array([50,50,50,50])
dfNi04 = "DF240317_Ni04_18cm_30MeV_calibration.Spe"


# channel_array_dgNi05 = np.array([1288, 2043]) #For få topper til å bruke i kalibrering
# energy_array_dgNi05 = np.array([510.9989461, 810.7593])
# energy_unc_array_dgNi05 = np.array([0.0000031, 0.002])
# dgNi05 = "DG240317_Ni05_18cm_30MeV_calibration.Spe"



# Defining which spectrum I want to calibrate
channel_array = channel_array_dfNi04 
energy_array = energy_array_dfNi04
energy_unc_array = energy_unc_array_dfNi04
subset_params = subset_params_dfNi04
spectrumname = dfNi04

# Load and plot the SPE file
filename = spectrumname
spectrum = []
with open(filename, 'r') as file:
    for line in file:
        spectrum.append(float(line))
counts = np.array(spectrum)
channels = np.arange(0, counts.size, 1)

plt.plot(channels, counts)
plt.title(f'{spectrumname}')
# plt.show()



# Define the Gaussian function
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Define the polynomial background function (quadratic)
def polynomial(x, a, b, c):
    return a * x**2 + b * x + c

# Combine the Gaussian and polynomial functions
def combined_function(x, A, mu, sigma, a, b, c):
    return gaussian(x, A, mu, sigma) + polynomial(x, a, b, c)



# Define initial guesses for fitting parameters (A, mu, sigma, a, b, c)
initial_guesses = []
for i in range(len(channel_array)):
    initial_guesses.append((100, channel_array[i], 10, 0, 1, 0))


# Perform combined Gaussian and polynomial fits
fits = []
channel_nr_fit_list = []
channel_nr_unc_list = []

for guess, subset in zip(initial_guesses, subset_params):
    # channels_subset = channels[(guess[1]-50):(guess[1]+50)]
    # counts_subset = counts[(guess[1]-50):(guess[1]+50)]
    channels_subset = channels[(guess[1]-subset):(guess[1]+subset)]
    counts_subset = counts[(guess[1]-subset):(guess[1]+subset)]
    # plt.plot(channels_subset, counts_subset, color = 'orchid')
    # plt.show()
    popt, cov = curve_fit(combined_function, channels_subset, counts_subset, p0=guess)
    fits.append(popt)
    print(popt[2])
    channel_nr_fit_list.append(popt[1])
    channel_nr_unc_list.append(popt[2]**2)

    # Calculate the residuals
    residuals = counts_subset - combined_function(channels_subset, *popt)

    # Calculate the degrees of freedom
    n_data = len(channels_subset)
    n_parameters = len(popt)
    degrees_of_freedom = n_data - n_parameters

    # Calculate the chi-squared value
    # chi_squared = np.sum((residuals / combined_function(channels_subset, *popt))**2) _________________________________________________XXXXXXXXXX____________________________
    chi_squared = np.sum(residuals**2 / combined_function(channels_subset, *popt))

    # Calculate the reduced chi-squared value
    reduced_chi_squared = chi_squared / degrees_of_freedom

    plt.plot(channels_subset, combined_function(channels_subset, *popt), label = f'reduced chi2: {reduced_chi_squared}')
plt.legend()
plt.show()




#Quadratic calibration 
polyline = np.linspace(0,7000, 700000)
calib_fit, cov_matrix = np.polyfit(channel_nr_fit_list, energy_array, 2, cov = True)
model = np.poly1d(calib_fit)
print("The energy calibration paramteters:")
print(model)
# Extract the standard deviations (uncertainties) of the fit parameters
calib_unc = np.sqrt(np.diag(cov_matrix))
print(f"The uncertainties of the energy calibration is {calib_unc}")


# Finding the chi2 of the calibration:
def chi_squared(observed, expected):
    return np.sum( (observed - expected)**2/expected )

chi_squared_calibration_fit = chi_squared(energy_array, model(channel_nr_fit_list))
print("\n")
print(f'Engcal:\n {calib_fit[2]}, \n {calib_fit[1]}, \n {calib_fit[0]}')

#Visualizing the calibration
plt.errorbar(channel_nr_fit_list, energy_array, xerr=channel_nr_unc_list, yerr = energy_unc_array, fmt='o', label = 'datapoints', color = 'skyblue')
plt.plot(polyline, model(polyline), color = 'hotpink', label = 'fit')
plt.legend()
plt.xlabel("Channel number")
plt.ylabel("Energy")
plt.text(3000, 10, f'Engcal:\n {calib_fit[2]} +- {calib_unc[2]:.2e} \n {calib_fit[1]} +- {calib_unc[1]:.2e}\n {calib_fit[0]} +- {calib_unc[0]:.2e} \n\n Chi2 = {chi_squared_calibration_fit}', fontsize=8, color='black', backgroundcolor='hotpink')
plt.title(f'{spectrumname}')
plt.show()








