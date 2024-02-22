import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd



dp = 0.977
# foil_name = 'Ni01'
# foil_list = ['Ni01', 'Zr01', 'Ti01', 'Ni02', 'Zr02', 'Ti02', 'Ni03', 'Zr03', 'Ti03', 'Ni04', 'Zr04', 'Ti04', 'Ni05', 'Zr05', 'Ti05']
# foil_list = ['Ti05']
# foil_list = ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05']
foil_list = ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05']
# foil_list = ['Ti01', 'Ti02', 'Ti03', 'Ti04', 'Ti05']



# colors = ['hotpink', ]
# Getting data from files_________________________________________________________________________________________________________________________________________
flux_file = f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.3f}_fluxes.csv'
csv_flux_data = pd.read_csv(flux_file)
# stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.3f}.csv')
foil_energy_data = {}

for foil_name in foil_list:
    # mu_energy = float(stack_df[stack_df['name']==foil_name]['mu_E'])

    target_flux_data = csv_flux_data.loc[csv_flux_data['name'] == foil_name]
    energy = target_flux_data.loc[:,'energy']
    flux = target_flux_data.loc[:,'flux']
    energy_from_stack_calc_unfiltered = np.array(energy.values.tolist())
    flux_from_stack_calc_unfiltered = np.array(flux.values.tolist()) #[MeV]
    energy_list= energy_from_stack_calc_unfiltered[energy_from_stack_calc_unfiltered >= 0.275]
    flux_list = flux_from_stack_calc_unfiltered[energy_from_stack_calc_unfiltered >= 0.275]
    mu_energy = np.sum(energy_list * flux_list)/np.sum(flux_list)


    # Calculating FWHM________________________________________________________________________________________________________________________________________________
    max_flux = np.max(flux_list)
    max_energy = energy_list[np.argmax(flux_list)]

    # Calculate the half-maximum flux
    half_max_flux = max_flux / 2.0

    # Find the indices where the flux crosses the half-maximum
    index_left = np.argmin(np.abs(flux_list[:np.argmax(flux_list)] - half_max_flux))
    index_right = np.argmin(np.abs(flux_list[np.argmax(flux_list):] - half_max_flux)) + np.argmax(flux_list)

    # Calculate the energies corresponding to these indices
    energy_left = energy_list[index_left]
    energy_right = energy_list[index_right]

    # Calculate the FWHM
    fwhm = energy_right - energy_left

    energy_min_unc = mu_energy - energy_list[index_left]
    energy_plus_unc = energy_list[index_right]-mu_energy

    # print(f'{foil_name}, fwhm = {fwhm:.2f}, E = {mu_energy:.2f} + {energy_plus_unc:.2f} - {energy_min_unc:.2f}')
    foil_energy_data[foil_name] = {'energy': mu_energy, 'min_unc': energy_min_unc, 'plus_unc': energy_plus_unc}

    # Plotting________________________________________________________________________________________________________________________________________________________
    # plt.plot(energy_list, flux_list, label = f'{foil_name}')
    plt.plot(energy_list, flux_list, label = f'{foil_name}, fwhm = {fwhm:.2f}')

    # plt.plot([max_energy, max_energy], [0, max_flux], color='grey', linestyle='--')
    plt.plot([mu_energy, mu_energy], [0, max_flux], color='grey', linestyle='--')
    plt.plot([energy_left, energy_right], [half_max_flux, half_max_flux] , color = 'black')
plt.xlabel('Energy (MeV)')
plt.ylabel('Reltive deutron flux')
plt.legend()
plt.show()

print(foil_energy_data)


