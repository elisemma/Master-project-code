import numpy as np  
import matplotlib.pyplot as plt 
import pandas as pd 
from foil_class_for_xs_calc import Foil 
from scipy.interpolate import splrep, splev
import os 
import csv
from x4i3 import exfor_manager, exfor_entry
from urllib.request import urlopen




def get_exfor_data(target, reaction, product): #f.exs: ('ZR-0', 'D,*', 'NB-90')
    # Pull data from EXFOR
    db = exfor_manager.X4DBManagerDefault()

    x = db.retrieve(target=target,reaction=reaction,quantity='SIG' )
    # Wildcards in 'reaction' seem to produce a TON of false positives, even when querying with 'product'
    # x = db.retrieve(target='BI-209',reaction='D,*',product='PO-209',quantity='SIG' )
    # query() just searches for data matching search parameters

    # Hold extracted data for plotting
    plot_Dict = {}
    for key in x.keys():
        entry = x[key]
        datasets = entry.getDataSets()

        # Make sure only one subentry is retrieved per entry!
        num_of_sub_subentries = len(list(datasets.keys()))

        if num_of_sub_subentries == 1:
            # We're all good...
            subentry = next(iter(datasets.values()))

        else:
            # Poorly-formatted EXFOR - multiple subentries for one entry
            print('Number of datasets found in entry', next(iter(datasets))[1][0:5], ': ', num_of_sub_subentries)

            for i in datasets.values():
                # Only select subentries leading to the specified product
                if product in str(i.reaction):
                    subentry = i
    

        # Pull metadata and stash into dictionary
        # print(subentry.getSimplified())
        # print(subentry.getSimplified().data)
        author_name = subentry.author[0].split('.',-1)[-1]
        year = subentry.year
        # print(author_name)
        # print(year)
        # print(type(subentry))
        plot_Dict[author_name] = (author_name, year, np.array(subentry.getSimplified().data, dtype=float), subentry.subent)

    return plot_Dict



def calc_xs(foil_name, reaction_product):
    #(foil_name, reaction_product, A0, A0_unc)
    A0_by_curie_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_curie.csv')
    if os.path.exists(f'./Calculated_A0/{foil_name}_A0_by_hand.csv'):
        A0_by_hand_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_hand.csv')
        A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)
    else:
        A0_concat_df = A0_by_curie_df

    A0_filtered_df = A0_concat_df[A0_concat_df['Isotope']==reaction_product]
    A0 = A0_filtered_df['A0'].iloc[0]
    A0_unc = A0_filtered_df['A0_unc'].iloc[0]

    foil = Foil(foil_name, reaction_product, A0, A0_unc)
    foil.assign_areal_dens_w_unc_percent()
    foil.assign_molar_mass()
    foil.calculate_decay_constant_w_unc()
    foil.assign_beam_current_w_unc()
    foil.calculate_xs_w_unc()

    xs = foil.calc_xs
    xs_unc = foil.calc_xs_unc
    return xs, xs_unc #[mb]



def get_talys_data(filename, foil_name):
    E_talys, xs_talys = np.loadtxt(f'./talys_calculations/talys_{foil_name[0:2]}/talys/'+filename, unpack = True, skiprows = 24)
    E_list = [0]+list(E_talys)
    xs_list = [0]+list(xs_talys)
    return E_list, xs_list




def plot_xs(reaction_product, Z, A, foil_list):

    energy_data = {'Ni01': {'energy': 27.337003000000003, 'min_unc': 0.6070030000000024, 'plus_unc': 0.6529969999999956}, 
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
                   'Ni05': {'energy': 2.2266116583208797, 'min_unc': 1.8966116583208799, 'plus_unc': 1.4033883416791202}, 
                   'Zr05': {'energy': 1.5675066930508672, 'min_unc': 1.2375066930508674, 'plus_unc': 0.7424933069491328}, 
                   'Ti05': {'energy': 0.7071428571428571, 'min_unc': 0.37714285714285717, 'plus_unc': 0.2828571428571429}}

    energy_list = [energy_data[foil]['energy'] for foil in foil_list]
    energy_min_unc_list = [energy_data[foil]['min_unc'] for foil in foil_list]
    energy_plus_unc_list = [energy_data[foil]['plus_unc'] for foil in foil_list]

    xs_list = []
    xs_unc_list = []

    for foil in foil_list:
        xs, xs_unc = calc_xs(foil, reaction_product)
        xs_list.append(xs)
        xs_unc_list.append(xs_unc)


    #Write calculated xs and energy w unc to file
    csv_file_path = f'./Calculated_xs/{reaction_product}_xs.csv'
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Foil', 'Energy', 'Energy_min_unc', 'Energy_plus_unc', 'xs', 'xs_unc'])
        for i in range(len(foil_list)):
            csv_writer.writerow([f'{foil_list[i]}', f'{energy_list[i]}', f'{energy_min_unc_list[i]}', f'{energy_plus_unc_list[i]}', f'{xs_list[i]}', f'{xs_unc_list[i]}'])

    #Get data from talys
    talys_file = f'rp0{Z}0{A}.tot'
    E_talys, xs_talys = get_talys_data(talys_file, foil_list[0])
    talys_spline = splrep(E_talys, xs_talys)
    energy_array = np.linspace(0,30,300)

    #Get data from exfor
    product = reaction_product[2:]+'-'+reaction_product[0:2]

    if foil_list[0][0:2]=='Ti':
        target = 'TI-0'
    elif foil_list[0][0:2]=='Zr':
        target = 'ZR-0'
    elif foil_list[0][0:2]=='Ni':
        target = 'NI-0'

    exfor_dict = get_exfor_data(target, 'D,*', product)


    #Plot results
    markers = ['.', '*', 'v', '^', '+', '<', '>', 's', 'h']
    grey_colors = ['dimgrey', 'darkgrey', 'lightgrey', 'silver', 'k', 'dimgrey', 'darkgrey', 'lightgrey', 'silver']
    k=0
    for index in exfor_dict:
        # Need to set up list of marker sizes to iterate over with k
        # print(k)
        if exfor_dict[index][2].shape[1] == 4:
            plt.errorbar(exfor_dict[index][2][:,0],1E3*exfor_dict[index][2][:,1], xerr=exfor_dict[index][2][:,2], yerr=1E3*exfor_dict[index][2][:,3], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])

        elif exfor_dict[index][2].shape[1] == 3:
            plt.errorbar(exfor_dict[index][2][:,0],1E3*exfor_dict[index][2][:,1], yerr=1E3*exfor_dict[index][2][:,2], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])

        elif exfor_dict[index][2].shape[1] >= 5:
            print('WARNING: Plotting',str(exfor_dict[index][2].shape[1])+'-column EXFOR data retrieved for subentry', exfor_dict[index][3]+', please make sure data look reasonable - column formatting is inconsistent for >4 columns.')
            plt.errorbar(exfor_dict[index][2][:,4],exfor_dict[index][2][:,5], xerr=exfor_dict[index][2][:,0], yerr=exfor_dict[index][2][:,6], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])
        k=k+1

    plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='gold', label='TALYS')
    plt.plot(energy_array, splev(energy_array, talys_spline)*0.75, linewidth=1, color='gold', label='TEST')
    plt.plot(energy_array, splev(energy_array, talys_spline)*1.2, linewidth=1, color='skyblue', label='Test')
    plt.plot(energy_array, splev(energy_array, talys_spline)*0.9, linewidth=1, color='mediumseagreen', label='test')
    plt.plot(energy_array, splev(energy_array, talys_spline)*1.1, linewidth=1, color='gold', label='test.')


    plt.errorbar(energy_list[0:len(foil_list)], xs_list, xerr=[energy_min_unc_list[0:len(foil_list)], energy_plus_unc_list[0:len(foil_list)]], yerr=xs_unc_list, marker='D', markersize=5, linewidth=2, linestyle='', color='deeppink', label='This work')
    plt.xlabel('Beam energy (MeV)')
    plt.ylabel('Cross section (mb)')
    plt.title(reaction_product)
    plt.legend()
    plt.xlim(0,30)
    plt.show()




#_________________________________Running the code_______________________________________

# plot_xs('96NB', 41, 96, ['Zr01', 'Zr02', 'Zr03', 'Zr04'])
plot_xs('90NB', 41, 90, ['Zr01', 'Zr02', 'Zr03', 'Zr04'])
# plot_xs('47SC', 21, 47, ['Ti01', 'Ti02', 'Ti03'])
# plot_xs('46SC', 21, 46, ['Ti01', 'Ti02', 'Ti03', 'Ti04'])
# plot_xs('48V', 23, 48, ['Ti01', 'Ti02', 'Ti03', 'Ti04'])
# plot_xs('52MN', 25, 52, ['Ni01', 'Ni02'])
# plot_xs('65NI', 28, 65, ['Ni01', 'Ni02', 'Ni03', 'Ni04'])
# plot_xs('57CO', 27, 57, ['Ni01', 'Ni02', 'Ni03'])
# plot_xs('60CU', 29, 60, ['Ni02', 'Ni03'])








