import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 



def calibration(channels, calibration_params):
    energies = calibration_params[0]*channels**2 + calibration_params[1]*channels + calibration_params[2]
    return energies
    



file = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/fit_peak_files/Ni_foils_analysis/peaks_fitted_manually/BE130217_Ni05_18cm_30MeV_peak_data_NO_PEAK_FITTED.Spe'

calibration_params_BE_Ni05 = [-0.26926665188089427, 0.3844465319472371, -1.2638762490788093e-07]


# Load and plot the SPE file
filename = file
spectrum = []
with open(filename, 'r') as file:
    for line in file:
        spectrum.append(float(line))
counts = np.array(spectrum)
channels = np.arange(0, counts.size, 1)


energies = calibration(channels, calibration_params_BE_Ni05 )

plt.plot(channels, counts)
plt.show()


 # Create histogram
# num_bins = 50
# histogram_values, bin_edges = np.histogram(counts, bins=num_bins, density=False)

plt.bar(channels, counts, color='hotpink', alpha=0.5, label='Data', width=channels[1]-channels[0])

plt.legend()
plt.title('Gaussian Peak with Linear Background')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()