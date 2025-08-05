import numpy as np 
import matplotlib.pyplot as plt 



xs_dict = {'Zr06': {'energy': 47.64,'xs_86': 39.1, 'xs_unc_86': 4.7, 'xs_87': 39.4, 'xs_unc_87': 2.9, 'xs_88': 353, 'xs_unc_88': 21, 'xs_90m': 6.20, 'xs_unc_90m': 0.59},
           'Zr07': {'energy': 41.1,'xs_86': 46.0, 'xs_unc_86': 4.8, 'xs_87': 26.9, 'xs_unc_87': 2.1, 'xs_88': 226, 'xs_unc_88': 19, 'xs_90m': 5.43, 'xs_unc_90m': 0.56},
           'Zr08': {'energy': 36.4,'xs_86': 33.1, 'xs_unc_86': 3.7, 'xs_87': 24.4, 'xs_unc_87': 3.1, 'xs_88': 80, 'xs_unc_88': 17, 'xs_90m': 5.27, 'xs_unc_90m': 0.74},
           'Zr09': {'energy': 31.0,'xs_86': 13.0, 'xs_unc_86': 1.8, 'xs_87': 31.8, 'xs_unc_87': 3.4, 'xs_88': 19.8, 'xs_unc_88': 1.2}}


# Calculate the ratios with uncertainties
for foil, values in xs_dict.items():
    # Calculate the ratio of 86/87
    ratio_86_87 = values['xs_86'] / values['xs_87']
    # Calculate the uncertainty of the ratio of 86/87 using propagation of uncertainties formula
    ratio_86_87_unc = ratio_86_87 * ((values['xs_unc_86'] / values['xs_86'])**2 + (values['xs_unc_87'] / values['xs_87'])**2)**0.5
    
    # Calculate the sum of cross sections for all other isotopes (88, 90m, etc.)
    all_yttrium_isotopes = sum(values[f'xs_{iso}'] for iso in ['87', '88', '90m'] if f'xs_{iso}' in values)
    
    # Calculate the ratio of 86 divided by the sum of all other isotopes
    ratio_86_all_yttrium = values['xs_86'] / all_yttrium_isotopes
    # Calculate the uncertainty of the ratio using propagation of uncertainties formula
    all_yttrium_uncertainties = sum(values[f'xs_unc_{iso}']**2 for iso in ['88', '90m'] if f'xs_{iso}' in values)
    ratio_86_all_yttrium_unc = ratio_86_all_yttrium * ((values['xs_unc_86'] / values['xs_86'])**2 + (all_yttrium_uncertainties / all_yttrium_isotopes**2))**0.5

    
    # Store the calculated ratios and their uncertainties back into the dictionary
    xs_dict[foil]['86_87_ratio'] = (ratio_86_87)
    xs_dict[foil]['86_87_ratio_unc'] = (ratio_86_87_unc)
    xs_dict[foil]['86_all_yttrium_ratio'] = (ratio_86_all_yttrium)
    xs_dict[foil]['86_all_yttrium_ratio_unc'] = (ratio_86_all_yttrium_unc)


print(xs_dict)


# Extract data for plotting
energies = [values['energy'] for values in xs_dict.values()]
ratios_86_87 = [values['86_87_ratio'] for values in xs_dict.values()]
uncertainties_86_87 = [values['86_87_ratio_unc'] for values in xs_dict.values()]
ratios_86_all_yttrium = [values['86_all_yttrium_ratio'] for values in xs_dict.values()]
uncertainties_86_all_yttrium = [values['86_all_yttrium_ratio_unc'] for values in xs_dict.values()]

# Plotting
plt.figure(figsize=(9,6))
plt.errorbar(energies, ratios_86_87, yerr=uncertainties_86_87, fmt='>', color='hotpink', capsize=5, label= r'$^{86}$Y/$^{87}$Y')
plt.errorbar(energies, ratios_86_all_yttrium, yerr=uncertainties_86_all_yttrium, fmt='<', color='skyblue', capsize=5, label= r'$^{86}$Y/ all observed yttrium isotopes ')
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Ratio', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
# plt.savefig('./Figures/ratios_of_86Y_vs_other_yttriums.pdf', dpi=600)
# plt.grid(True)
plt.show()