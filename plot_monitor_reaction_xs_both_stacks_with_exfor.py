import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from foil_class import Foil
import os 
from scipy.interpolate import PchipInterpolator
import re



def get_IAEA_monitro_xs(target_material, reaction_product):
    filename = './Monitor_cross_section_data/IAEA_monitor_xs_' + target_material + '_dx_' + reaction_product + '.txt'
    E_mon_list = []
    xs_list = []
    xs_unc_list = []
    with open(filename) as file:
        lines = file.readlines()[6:-1]
        
        for line in lines:
            words = line.split()
            E_mon_list.append(float(words[0]))
            xs_list.append(float(words[1]))
            xs_unc_list.append(float(words[2]))
        file.close()
    interp_xs = PchipInterpolator(E_mon_list, xs_list)
    interp_unc_xs = PchipInterpolator(E_mon_list, xs_unc_list)

    # return E_mon_list, xs_list, xs_unc_list
    return interp_xs, interp_unc_xs



#_________________________Energy info from another file________________________
foil_energy_data = {'Ni01': {'energy': 27.337003000000003, 'min_unc': 0.6070030000000024, 'plus_unc': 0.6529969999999956}, 
                    'Zr01': {'energy': 26.375581480000005, 'min_unc': 0.605581480000005, 'plus_unc': 0.5944185199999943}, 
                    'Ti01': {'energy': 25.518248420000003, 'min_unc': 0.5882484200000029, 'plus_unc': 0.6117515799999964}, 
                    'Ni02': {'energy': 20.392268865, 'min_unc': 0.7422688649999998, 'plus_unc': 0.7577311350000002}, 
                    'Zr02': {'energy': 19.1745053, 'min_unc': 0.7845052999999993, 'plus_unc': 0.7754946999999994}, 
                    'Ti02': {'energy': 18.058134579999997, 'min_unc': 0.8081345799999973, 'plus_unc': 0.8118654200000037}, 
                    'Ni03': {'energy': 13.9680305, 'min_unc': 0.9780305000000009, 'plus_unc': 1.0619695}, 
                    'Zr03': {'energy': 12.309193775999999, 'min_unc': 1.059193775999999, 'plus_unc': 1.1008062240000012}, 
                    'Ti03': {'energy': 10.706294196000002, 'min_unc': 1.1362941960000015, 'plus_unc': 1.2037058039999984}, 
                    'Ni04': {'energy': 8.648899379944393, 'min_unc': 1.2988993799443929, 'plus_unc': 1.5211006200556056}, 
                    'Zr04': {'energy': 6.082350534817555, 'min_unc': 1.4923505348175548, 'plus_unc': 1.8676494651824447}, 
                    'Ti04': {'energy': 3.6998866038357203, 'min_unc': 2.2298866038357206, 'plus_unc': 2.090113396164279}, 
                    # 'Ni05': {'energy': 2.2266116583208797, 'min_unc': 1.8966116583208799, 'plus_unc': 1.4033883416791202}, 
                    # 'Zr05': {'energy': 1.5675066930508672, 'min_unc': 1.2375066930508674, 'plus_unc': 0.7424933069491328},
                    # 'Ti05': {'energy': 0.7071428571428571, 'min_unc': 0.37714285714285717, 'plus_unc': 0.2828571428571429},
                    'Fe01': {'energy': 48.25110370000002, 'min_unc': 0.9011037000000144, 'plus_unc': 0.8988962999999828}, 
                    'Zr06': {'energy': 47.66121943333334, 'min_unc': 0.9112194333333434, 'plus_unc': 0.8887805666666537}, 
                    'Ti06': {'energy': 47.13022093333334, 'min_unc': 0.8802209333333408, 'plus_unc': 0.9197790666666563}, 
                    'Fe02': {'energy': 41.876494533333336, 'min_unc': 1.0264945333333273, 'plus_unc': 1.073505466666667}, 
                    'Zr07': {'energy': 41.21671893333334, 'min_unc': 1.066718933333334, 'plus_unc': 1.0332810666666603}, 
                    'Ti08': {'energy': 40.62009246666667, 'min_unc': 1.0700924666666722, 'plus_unc': 1.0299075333333363}, 
                    'Fe03': {'energy': 37.24984646666667, 'min_unc': 1.099846466666662, 'plus_unc': 1.1001535333333408}, 
                    'Zr08': {'energy': 36.523388133333334, 'min_unc': 1.173388133333333, 'plus_unc': 1.1266118666666713}, 
                    'Ti09': {'energy': 35.86358033333333, 'min_unc': 1.1135803333333314, 'plus_unc': 1.1864196666666658}, 
                    'Fe04': {'energy': 32.115906233333334, 'min_unc': 1.2659062333333324, 'plus_unc': 1.2340937666666676}, 
                    'Zr09': {'energy': 31.29907423333333, 'min_unc': 1.249074233333328, 'plus_unc': 1.2509257666666684}, 
                    'Ti10': {'energy': 30.5522245, 'min_unc': 1.3022245000000012, 'plus_unc': 1.2977755000000002}, 
                    'Fe05': {'energy': 26.2370075, 'min_unc': 1.4870075000000007, 'plus_unc': 1.5129924999999993}, 
                    'Zr10': {'energy': 25.275796966666668, 'min_unc': 1.5257969666666682, 'plus_unc': 1.5742030333333332}, 
                    'Ti11': {'energy': 24.390982466666664, 'min_unc': 1.5409824666666623, 'plus_unc': 1.5590175333333391}}



data_by_reaction_p1 = {'Ni_dx_56CO': {'calc_xs': [np.float64(10.532442181808857), np.float64(8.408885213771274), np.float64(20.95950924976207), np.float64(31.352732387066354)], 
                                'calc_xs_unc': [np.float64(1.9771680973277075), np.float64(1.4363308990558326), np.float64(1.6295246024844563), np.float64(2.6836661486484665)], 
                                'energy': [np.float64(27.35133924), np.float64(20.450330565), np.float64(14.083908879999997), np.float64(8.848522071305181)]}, 
                       
                       'Ni_dx_58CO': {'calc_xs': [np.float64(137.87837417453406), np.float64(225.35754150734147), np.float64(176.8409492561857), np.float64(48.99062836248532)], 
                                 'calc_xs_unc': [np.float64(25.34306071111694), np.float64(50.56467095954041), np.float64(13.412564222612255), np.float64(18.45017123351097)], 
                                 'energy': [np.float64(27.35133924), np.float64(20.450330565), np.float64(14.083908879999997), np.float64(8.848522071305181)]}, 
                                 
                        'Ni_dx_61CU': {'calc_xs': [np.float64(15.42066162933613), np.float64(15.774345744614937), np.float64(25.155185910486416), np.float64(59.46452785548102)], 
                                 'calc_xs_unc': [np.float64(2.721863036809532), np.float64(2.708878545346285), np.float64(1.9531371664113437), np.float64(5.189430812285558)], 
                                 'energy': [np.float64(27.35133924), np.float64(20.450330565), np.float64(14.083908879999997), np.float64(8.848522071305181)]}, 
                                 
                        'Ti_dx_48V': {'calc_xs': [np.float64(53.28131917888242), np.float64(71.53510820381517), np.float64(86.95500752412562), np.float64(120.32118113088468), np.float64(211.09294415195407), np.float64(186.55464931381064), np.float64(325.8026945250252), np.float64(153.00128587401133), np.float64(17.443017872325555)], 
                                'calc_xs_unc': [np.float64(4.004285252699319), np.float64(2.215160347196869), np.float64(2.3840605321553276), np.float64(5.720374572911209), np.float64(2.7645225694942264), np.float64(7.584885271590837), np.float64(11.540288232419037), np.float64(6.703507158349622), np.float64(8.539100018271506)], 
                                'energy': [np.float64(47.09795026666667), np.float64(40.50592943333334), np.float64(35.68120613333333), np.float64(30.28116373333334), np.float64(23.99016156666667), np.float64(25.542917359999997), np.float64(18.134367400000002), np.float64(10.866410484000001), np.float64(3.9274296448168275)]}, 
                                
                        'Ti_dx_46SC': {'calc_xs': [np.float64(71.23392223363822), np.float64(47.64071433003223), np.float64(31.38176508759792), np.float64(25.349146020414114), np.float64(23.97580225065507), np.float64(23.856498720283238), np.float64(25.26628579683382), np.float64(33.68948658218453), np.float64(3.136104517890606)], 
                                 'calc_xs_unc': [np.float64(5.412575994612804), np.float64(1.4596181121256113), np.float64(0.8335479522980904), np.float64(1.1962400980743435), np.float64(0.3986503440442611), np.float64(0.9951421920697885), np.float64(0.8973800163182105), np.float64(1.4889929111037854), np.float64(1.537752776847068)], 
                                 'energy': [np.float64(47.09795026666667), np.float64(40.50592943333334), np.float64(35.68120613333333), np.float64(30.28116373333334), np.float64(23.99016156666667), np.float64(25.542917359999997), np.float64(18.134367400000002), np.float64(10.866410484000001), np.float64(3.9274296448168275)]},
                        
                        'Fe_dx_56CO': {'calc_xs': [np.float64(42.80199907072412), np.float64(53.53658424185413), np.float64(67.42712662562114), np.float64(95.89102219731262), np.float64(165.9191705527795)], 
                                 'calc_xs_unc': [np.float64(3.0335331140655164), np.float64(3.843651004670334), np.float64(4.898118323109198), np.float64(6.741815080623053), np.float64(11.665288091765959)], 
                                 'energy': [np.float64(48.2318781), np.float64(41.779047799999994), np.float64(37.08828046666665), np.float64(31.872687066666668), np.float64(25.879117466666663)]}}




#____________________________Getting monitor xs__________________________________
E_mon_array = np.linspace(0,50,100000)

interp_xs_Ti_dx_46SC, interp_unc_xs_Ti_dx_46SC = get_IAEA_monitro_xs('Ti','46SC')
interp_xs_Ti_dx_48V, interp_unc_xs_Ti_dx_48V = get_IAEA_monitro_xs('Ti','48V')
interp_xs_Ni_dx_56CO, interp_unc_xs_Ni_dx_56CO = get_IAEA_monitro_xs('Ni','56CO')
interp_xs_Ni_dx_58CO, interp_unc_xs_Ni_dx_58CO = get_IAEA_monitro_xs('Ni','58CO')
interp_xs_Ni_dx_61CU, interp_unc_xs_Ni_dx_61CU = get_IAEA_monitro_xs('Ni','61CU')
interp_xs_Fe_dx_56CO, interp_unc_xs_Fe_dx_56CO = get_IAEA_monitro_xs('Fe','56CO')




Ni_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Ni')]
Ni_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Ni')]
Ti_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Ti')]
Ti_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Ti')]
Fe_energy_min_unc_list = [foil_energy_data[foil]['min_unc'] for foil in foil_energy_data if foil.startswith('Fe')]
Fe_energy_plus_unc_list = [foil_energy_data[foil]['plus_unc'] for foil in foil_energy_data if foil.startswith('Fe')]


calc_xs_Ni_dx_56CO_p1 = data_by_reaction_p1['Ni_dx_56CO']['calc_xs']
calc_xs_unc_Ni_dx_56CO_p1 = data_by_reaction_p1['Ni_dx_56CO']['calc_xs_unc']
energies_Ni_dx_56CO_p1 = data_by_reaction_p1['Ni_dx_56CO']['energy']

calc_xs_Ni_dx_58CO_p1 = data_by_reaction_p1['Ni_dx_58CO']['calc_xs']
calc_xs_unc_Ni_dx_58CO_p1 = data_by_reaction_p1['Ni_dx_58CO']['calc_xs_unc']
energies_Ni_dx_58CO_p1 = data_by_reaction_p1['Ni_dx_58CO']['energy']

calc_xs_Ni_dx_61CU_p1 = data_by_reaction_p1['Ni_dx_61CU']['calc_xs']
calc_xs_unc_Ni_dx_61CU_p1 = data_by_reaction_p1['Ni_dx_61CU']['calc_xs_unc']
energies_Ni_dx_61CU_p1 = data_by_reaction_p1['Ni_dx_61CU']['energy']

calc_xs_Ti_dx_46SC_p1 = data_by_reaction_p1['Ti_dx_46SC']['calc_xs']
calc_xs_unc_Ti_dx_46SC_p1 = data_by_reaction_p1['Ti_dx_46SC']['calc_xs_unc']
energies_Ti_dx_46SC_p1 = data_by_reaction_p1['Ti_dx_46SC']['energy']

calc_xs_Ti_dx_48V_p1 = data_by_reaction_p1['Ti_dx_48V']['calc_xs']
calc_xs_unc_Ti_dx_48V_p1 = data_by_reaction_p1['Ti_dx_48V']['calc_xs_unc']
energies_Ti_dx_48V_p1 = data_by_reaction_p1['Ti_dx_48V']['energy']

calc_xs_Fe_dx_56CO_p1 = data_by_reaction_p1['Fe_dx_56CO']['calc_xs']
calc_xs_unc_Fe_dx_56CO_p1 = data_by_reaction_p1['Fe_dx_56CO']['calc_xs_unc'] 
energies_Fe_dx_56CO_p1 = data_by_reaction_p1['Fe_dx_56CO']['energy']




# ______________________________Plotting________________________________________

#Get data from exfor
markers = ['.', '*', 'v', '^', '+', '<', '>', 's', 'h',     '.', '*', 'v', '^', '+', '<', '>', 's', 'h']
grey_colors = ['dimgrey', 'darkgrey', 'lightgrey', 'silver', 'k', 'dimgrey', 'darkgrey', 'lightgrey', 'silver',     'silver', 'lightgrey', 'darkgrey', 'dimgrey', 'k','silver', 'lightgrey', 'dimgrey', 'darkgrey']



plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ni_dx_56CO/'
ref_numbers_Ni_dx_56Co = ['[45]', '[41]', '[93]', '[42]', '[94]', '[44]']
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year}) {ref_numbers_Ni_dx_56Co[i]}"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_56CO_p0, calc_xs_56CO_p0, xerr=[Ni_energy_min_unc_list[:-1], Ni_energy_plus_unc_list[:-1]], yerr=calc_xs_unc_56CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ni_dx_56CO(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ni_dx_56CO(E_mon_array)-interp_unc_xs_Ni_dx_56CO(E_mon_array), interp_xs_Ni_dx_56CO(E_mon_array)+interp_unc_xs_Ni_dx_56CO(E_mon_array), color='lavender')
plt.errorbar(energies_Ni_dx_56CO_p1, calc_xs_Ni_dx_56CO_p1, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_Ni_dx_56CO_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{56}$Co', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.ylim(0,60)
plt.savefig('./Figures/PhD/xs_mon_Ni_dx_56Co.pdf', dpi=600)
plt.show()



plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ni_dx_58CO/'
ref_numbers_Ni_dx_58Co = ['[41]', '[42]', '[64]', '[44]']
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year}) {ref_numbers_Ni_dx_58Co[i]}"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_58CO_p0, calc_xs_58CO_p0, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_58CO_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ni_dx_58CO(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ni_dx_58CO(E_mon_array)-interp_unc_xs_Ni_dx_58CO(E_mon_array), interp_xs_Ni_dx_58CO(E_mon_array)+interp_unc_xs_Ni_dx_58CO(E_mon_array), color='lavender')
plt.errorbar(energies_Ni_dx_58CO_p1, calc_xs_Ni_dx_58CO_p1, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_Ni_dx_58CO_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{58}$Co', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.ylim(0)
plt.savefig('./Figures/PhD/xs_mon_Ni_dx_58Co.pdf', dpi=600)
plt.show()



plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ni_dx_61CU/'
ref_numbers_Ni_dx_61Cu = ['[98]', '[41]', '[47]', '[48]', '[95]', '[44]', '[96]']
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year}) {ref_numbers_Ni_dx_61Cu[i]}"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_61CU_p0, calc_xs_61CU_p0, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_61CU_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ni_dx_61CU(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ni_dx_61CU(E_mon_array)-interp_unc_xs_Ni_dx_61CU(E_mon_array), interp_xs_Ni_dx_61CU(E_mon_array)+interp_unc_xs_Ni_dx_61CU(E_mon_array), color='lavender')
plt.errorbar(energies_Ni_dx_61CU_p1, calc_xs_Ni_dx_61CU_p1, xerr=[Ni_energy_min_unc_list, Ni_energy_plus_unc_list], yerr=calc_xs_unc_Ni_dx_61CU_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ni(d,x)$^{61}$Cu', fontsize=14)
plt.legend()
plt.xlim(0,30)
plt.ylim(0)
plt.savefig('./Figures/PhD/xs_mon_Ni_dx_61Cu.pdf', dpi=600)
plt.show()



plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ti_dx_46SC/'
ref_numbers_Ti_dx_46Sc = ['[50]', '[51]', '[52]', '[97]', '[54]', '[56]']
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year}) {ref_numbers_Ti_dx_46Sc[i]}"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_46SC_p0, calc_xs_46SC_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_46SC_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ti_dx_46SC(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ti_dx_46SC(E_mon_array)-interp_unc_xs_Ti_dx_46SC(E_mon_array), interp_xs_Ti_dx_46SC(E_mon_array)+interp_unc_xs_Ti_dx_46SC(E_mon_array), color='lavender')
plt.errorbar(energies_Ti_dx_46SC_p1, calc_xs_Ti_dx_46SC_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_Ti_dx_46SC_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ti(d,x)$^{46}$Sc', fontsize=14)
plt.legend()
plt.xlim(0,50)
plt.ylim(0)
plt.savefig('./Figures/PhD/xs_mon_Ti_dx_46Sc.pdf', dpi=600)
plt.show()



plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Ti_dx_48V/'
ref_numbers_Ti_dx_48V = ['[98]', '[50]', '[51]', '[52]']
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year}) {ref_numbers_Ti_dx_48V[i]}"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
# plt.errorbar(energies_48V_p0, calc_xs_48V_p0, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_48V_p0, marker='d', markersize=5, linestyle='', color='deepskyblue', label='p0')
plt.plot(E_mon_array, interp_xs_Ti_dx_48V(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Ti_dx_48V(E_mon_array)-interp_unc_xs_Ti_dx_48V(E_mon_array), interp_xs_Ti_dx_48V(E_mon_array)+interp_unc_xs_Ti_dx_48V(E_mon_array), color='lavender')
plt.errorbar(energies_Ti_dx_48V_p1, calc_xs_Ti_dx_48V_p1, xerr = [Ti_energy_min_unc_list, Ti_energy_plus_unc_list], yerr=calc_xs_unc_Ti_dx_48V_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Ti(d,x)$^{48}$V', fontsize=14)
plt.legend()
plt.xlim(0,50)
plt.ylim(0)
plt.savefig('./Figures/PhD/xs_mon_Ti_dx_48V.pdf', dpi=600)
plt.show()



plt.figure(figsize=(8, 6))
path = f'./Monitor_cross_section_data/Recommended_monitor_xs_Fe_dx_56CO/'
ref_numbers_Fe_dx_56Co = ['[66]', '[57]', '[98]', '[45]', '[63]', '[58]', '[59]']
 # List all files in the folder
files = os.listdir(path)
# Loop through each file
for i, file_name in enumerate(files):
    file_path = os.path.join(path, file_name)
    data = np.loadtxt(file_path, comments=["#", "//"])
    energy = data[:, 0]
    energy_unc = data[:, 1]
    cross_section = data[:, 2] * 1e3
    cross_section_unc = data[:, 3] * 1e3
    # Initialize variables to store author name and year
    author_name = ""
    year = ""
    # Read the file to find author name and year
    with open(file_path, 'r') as file:
        for line in file:
            if '# ' in line:  # Check if the line contains a comment
                # Use regular expression to find the author name and year
                match = re.search(r'(\d{4}),([^#]+)', line)
                if match:
                    year = match.group(1)  # Extract the year
                    author_name = match.group(2).strip()  # Extract the author name
                    # Remove unwanted single-letter followed by a dot and the plus sign
                    author_name = re.sub(r'\b[A-Z]\.\b|\+', '', author_name)
                    break  # Stop after finding the first occurrence of name and year
    # Create the label
    label = f"{author_name} ({year}) {ref_numbers_Fe_dx_56Co[i]}"
    plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)
plt.plot(E_mon_array, interp_xs_Fe_dx_56CO(E_mon_array), color='cornflowerblue', label='IAEA')
plt.fill_between(E_mon_array, interp_xs_Fe_dx_56CO(E_mon_array)-interp_unc_xs_Fe_dx_56CO(E_mon_array), interp_xs_Fe_dx_56CO(E_mon_array)+interp_unc_xs_Fe_dx_56CO(E_mon_array), color='lavender')
plt.errorbar(energies_Fe_dx_56CO_p1, calc_xs_Fe_dx_56CO_p1, xerr=[Fe_energy_min_unc_list, Fe_energy_plus_unc_list], yerr=calc_xs_unc_Fe_dx_56CO_p1, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
plt.ylabel('Cross Section (mb)', fontsize=14)
plt.title(r'$^{nat}$Fe(d,x)$^{56}$Co', fontsize=14)
plt.legend()
plt.xlim(0,50)
plt.ylim(0)
plt.savefig('./Figures/PhD/xs_mon_Fe_dx_56Co.pdf', dpi=600)
plt.show()




