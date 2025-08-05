import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy import odr


# Define the Gaussian function
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Define the polynomial background function (quadratic)
def polynomial(x, a, b, c):
    return a * x**2 + b * x + c

# Combine the Gaussian and polynomial functions
def combined_function(x, A, mu, sigma, a, b, c):
    return gaussian(x, A, mu, sigma) + polynomial(x, a, b, c)


path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/fit_peak_files/Ni_foils_analysis/peaks_fitted_manually/'

channel_array_daNi01 = np.array([1771]) #Using these values to make a calibration for all the Ni foils (DF and DG)
energy_array_daNi01 = np.array([1771])
energy_unc_array_daNi01 = np.array([0.002])
subset_params_daNi01 = np.array([50])
daNi01 = "DA200317_Ni01_18cm_30MeV_peak_data_NO_PEAK_FITTED.Spe"


channel_array_dbNi02 = np.array([1771]) #Using these values to make a calibration for all the Ni foils (DF and DG)
energy_array_dbNi02 = np.array([1771])
energy_unc_array_dbNi02 = np.array([0.002])
subset_params_dbNi02 = np.array([50])
dbNi02 = "DB200317_Ni02_18cm_30MeV_peak_data_NO_PEAK_FITTED.Spe"


channel_array_beNi05 = np.array([1332, 2114]) #Using these values to make a calibration for all the Ni foils (DF and DG)
energy_array_beNi05 = np.array([511, 810])
energy_unc_array_beNi05 = np.array([1, 1])
subset_params_beNi05 = np.array([50,50])
beNi05 = "BE130217_Ni05_18cm_30MeV_peak_data_NO_PEAK_FITTED.Spe"


# Defining which spectrum I want to calibrate
channel_array = channel_array_beNi05 
energy_array = energy_array_beNi05
energy_unc_array = energy_unc_array_beNi05
subset_params = subset_params_beNi05
spectrumname = beNi05

# Load and plot the SPE file
filename = path + spectrumname
spectrum = []
with open(filename, 'r') as file:
    for line in file:
        spectrum.append(float(line))
counts = np.array(spectrum)
channels = np.arange(0, counts.size, 1)

plt.plot(channels, counts)
plt.title(f'{spectrumname}')
# plt.show()



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
    channel_nr_unc_list.append(np.sqrt(cov[1,1]))

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

    plt.plot(channels_subset, combined_function(channels_subset, *popt), label = f'reduced chi2: {reduced_chi_squared}', linestyle = 'dashed')
plt.legend()
plt.show()




#Quadratic calibration 

def quad_poly(B, x):
    return B[0]*x**2 + B[1]*x + B[2]

data = odr.RealData(channel_nr_fit_list,energy_array,sx = np.array(channel_nr_unc_list), sy = energy_unc_array)
calib_model = odr.Model(quad_poly)
odr = odr.ODR(data, calib_model,beta0 =[0,1,0])
out = odr.run()#Run the regression
calib_params = out.beta#out.beta is a list which contains the values of the parameters popped out from the fitting
calib_params_unc = out.sd_beta#out.sd_beta is a list which contains the values of the errors of the parameters
print("The energy calibration paramteters:")
print(calib_params)
# Extract the standard deviations (uncertainties) of the fit parameters
print(f"The uncertainties of the energy calibration is {calib_params_unc}")

polyline = np.linspace(0,7000, 700000)
model = np.poly1d(calib_params)

print("\n")
print(f'Engcal:\n {calib_params[2]}, \n {calib_params[1]}, \n {calib_params[0]}')




#Visualizing the calibration
plt.errorbar(channel_nr_fit_list, energy_array, xerr=channel_nr_unc_list, yerr = energy_unc_array, fmt='o', label = 'datapoints', color = 'skyblue')
plt.plot(polyline, model(polyline), color = 'hotpink', label = 'fit')
plt.legend()
plt.xlabel("Channel number")
plt.ylabel("Energy")
plt.text(3000, 10, f'Engcal:\n {calib_params[2]} +- {calib_params_unc[2]:.2e} \n {calib_params[1]} +- {calib_params_unc[1]:.2e}\n {calib_params[0]} +- {calib_params_unc[0]:.2e} \n\n Reduced Chi2 = {out.res_var}', fontsize=8, color='black', backgroundcolor='hotpink')
plt.title(f'{spectrumname}')
plt.show()