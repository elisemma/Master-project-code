import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd





# foil_list_30MeV = ['Ni01', 'Zr01', 'Ti01', 'Ni02', 'Zr02', 'Ti02', 'Ni03', 'Zr03', 'Ti03', 'Ni04', 'Zr04', 'Ti04', 'Ni05', 'Zr05', 'Ti05']
# foil_list_50MeV = ['Fe01', 'Zr06', 'Ti06', 'Fe02', 'Zr07', 'Ti08', 'Fe03', 'Zr08', 'Ti09', 'Fe04', 'Zr09', 'Ti10', 'Fe05', 'Zr10', 'Ti11']

foil_list_30MeV = ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05']
foil_list_50MeV = ['Zr06','Zr07', 'Zr08', 'Zr09', 'Zr10']

colors = ['mediumseagreen', 'cornflowerblue', 'royalblue', 'm', 'violet', 'hotpink', 'sandybrown', 'gold', 'limegreen', 'deepskyblue']
i=0

# Getting data from files_________________________________________________________________________________________________________________________________________
foil_energy_data = {}
plt.figure(figsize=(10, 7))

for foil_name in foil_list_30MeV:
    # mu_energy = float(stack_df[stack_df['name']==foil_name]['mu_E'])
    dp = 0.972
    flux_file = f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.3f}_fluxes.csv'
    csv_flux_data = pd.read_csv(flux_file)

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
    plt.plot(energy_list, flux_list, color = colors[i], label = f'{foil_name}, fwhm = {fwhm:.2f}')

    # plt.plot([max_energy, max_energy], [0, max_flux], color='grey', linestyle='--')
    plt.plot([mu_energy, mu_energy], [0, max_flux], color='grey', linestyle='--')
    plt.plot([energy_left, energy_right], [half_max_flux, half_max_flux] , color = 'black')
    i+=1



for foil_name in foil_list_50MeV:
    dp = 0.990
    flux_file = f'./50MeV_stack_analysis/Stack_calculations/stack_50MeV_dp_{dp:.3f}_fluxes.csv'
    csv_flux_data = pd.read_csv(flux_file)
    # mu_energy = float(stack_df[stack_df['name']==foil_name]['mu_E'])

    target_flux_data = csv_flux_data.loc[csv_flux_data['name'] == foil_name]
    energy = target_flux_data.loc[:,'energy']
    flux = target_flux_data.loc[:,'flux']
    energy_list = np.array(energy.values.tolist())
    flux_list = np.array(flux.values.tolist()) #[MeV]
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

    foil_energy_data[foil_name] = {'energy': mu_energy, 'min_unc': energy_min_unc, 'plus_unc': energy_plus_unc}

    # Plotting________________________________________________________________________________________________________________________________________________________
    # plt.plot(energy_list, flux_list, label = f'{foil_name}')
    plt.plot(energy_list, flux_list, color = colors[i], label = f'{foil_name}, fwhm = {fwhm:.2f}')

    # plt.plot([max_energy, max_energy], [0, max_flux], color='grey', linestyle='--')
    plt.plot([mu_energy, mu_energy], [0, max_flux], color='grey', linestyle='--')
    plt.plot([energy_left, energy_right], [half_max_flux, half_max_flux] , color = 'black')
    i += 1







plt.xlabel('Deuteron Energy (MeV)', fontsize=12)
plt.ylabel('Reltive Deutron Flux', fontsize=12)
plt.title('Energy distribution for Zr foils', fontsize=12)
plt.legend(fontsize=11)
# plt.savefig('./Figures/Zr_energies.pdf', dpi=600)
plt.show()



