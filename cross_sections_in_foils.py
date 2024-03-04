import numpy as np  
import matplotlib.pyplot as plt 
import pandas as pd 
from foil_class_for_xs_calc import Foil 
from scipy.interpolate import splrep, splev
import os 
import csv
from x4i3 import exfor_manager, exfor_entry
from urllib.request import urlopen
from scipy.interpolate import PchipInterpolator





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



def calc_xs_from_A0(foil_name, reaction_product): #fungerer ikke med denne klassen
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



def calc_xs_from_R(foil_name, reaction_product): #fungerer ikke med denne klassen
    #(foil_name, reaction_product, R, R_unc)
    R_by_curie_df = pd.read_csv(f'./Calculated_R/{foil_name}_R_by_curie.csv')
    print(R_by_curie_df)
    R_filtered_df = R_by_curie_df[R_by_curie_df['Isotope']==reaction_product]
    R = R_filtered_df['R'].iloc[0]
    R_unc = R_filtered_df['R_unc'].iloc[0]

    foil = Foil(foil_name, reaction_product, R, R_unc)
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


def get_alice_data(target_material, Z, A):
    # Initialize lists to store energy and cross-section data
    file_path = f'./alice_calculations/plot_nat{target_material}_dx'
    energies = []
    cross_sections = []

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        # Skip lines until reaching the cross-section data
        for line in file:
            if 'Ebeam' in line:
                break
        # Read cross-section data for the specified isotope
        for line in file:
            if line.strip():  # Check if line is not empty
                data = line.split()
                try:
                    if int(float(data[1])) == Z and int(float(data[2])) == A:
                        energies.append(float(data[0]))
                        cross_sections.append(float(data[3]))
                except ValueError:
                    pass  # Skip lines that can't be converted to integers
    
    return energies, cross_sections




def plot_xs(reaction_product, Z, A, foil_list, write_csv=False, save_fig=False, use_R_to_calc_xs=True):

    energy_data = {'Ni01': {'energy': 27.35133924, 'min_unc': 0.6213392400000011, 'plus_unc': 0.6386607599999969}, 
                   'Zr01': {'energy': 26.395312639999993, 'min_unc': 0.6253126399999935, 'plus_unc': 0.634687360000008}, 
                   'Ti01': {'energy': 25.542917359999997, 'min_unc': 0.6129173599999973, 'plus_unc': 0.6470826400000043}, 
                   'Ni02': {'energy': 20.450330565, 'min_unc': 0.7403305650000007, 'plus_unc': 0.7596694349999993}, 
                   'Zr02': {'energy': 19.241718260000003, 'min_unc': 0.7917182600000032, 'plus_unc': 0.768281739999999}, 
                   'Ti02': {'energy': 18.134367400000002, 'min_unc': 0.7643674000000047, 'plus_unc': 0.7956325999999976}, 
                   'Ni03': {'energy': 14.083908879999997, 'min_unc': 0.973908879999998, 'plus_unc': 1.0060911200000024}, 
                   'Zr03': {'energy': 12.445472363999999, 'min_unc': 1.015472363999999, 'plus_unc': 1.0845276360000007}, 
                   'Ti03': {'energy': 10.866410484000001, 'min_unc': 1.116410484000001, 'plus_unc': 1.2235895159999988}, 
                   'Ni04': {'energy': 8.848522071305181, 'min_unc': 1.318522071305182, 'plus_unc': 1.4414779286948178}, 
                   'Zr04': {'energy': 6.346733794120523, 'min_unc': 1.456733794120523, 'plus_unc': 1.8432662058794769}, 
                   'Ti04': {'energy': 3.9274296448168275, 'min_unc': 2.1574296448168275, 'plus_unc': 2.1625703551831723}, 
                   'Ni05': {'energy': 2.3286219274448703, 'min_unc': 1.9986219274448704, 'plus_unc': 1.4813780725551298}, 
                   'Zr05': {'energy': 1.6333420201661875, 'min_unc': 1.3033420201661876, 'plus_unc': 0.8566579798338128}, 
                   'Ti05': {'energy': 1.3297790055248617, 'min_unc': 0.9997790055248619, 'plus_unc': -0.3997790055248617}}

    energy_list = [energy_data[foil]['energy'] for foil in foil_list]
    energy_min_unc_list = [energy_data[foil]['min_unc'] for foil in foil_list]
    energy_plus_unc_list = [energy_data[foil]['plus_unc'] for foil in foil_list]

    xs_list = []
    xs_unc_list = []

    for foil in foil_list:
        if use_R_to_calc_xs==True:
            xs, xs_unc = calc_xs_from_R(foil, reaction_product)
            xs_list.append(xs)
            xs_unc_list.append(xs_unc)
        else:
            xs, xs_unc = calc_xs_from_A0(foil, reaction_product)
            xs_list.append(xs)
            xs_unc_list.append(xs_unc)


    #Write calculated xs and energy w unc to file
    if write_csv == True:
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


    #Get data from alice
    target_material = foil_list[0][0:2]
    E_alice, xs_alice = get_alice_data(target_material, Z, A)
    E_alice = [0]+E_alice
    xs_alice = [0]+xs_alice

    # alice_spline = splrep(E_alice, xs_alice)
    alice_spline = PchipInterpolator(E_alice, xs_alice)

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

    # plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='gold', label='TALYS')
    # plt.plot(energy_array, splev(energy_array, alice_spline), linewidth=1, color='plum', linestyle='--', label='ALICE')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.2, linewidth=1, color='skyblue', label='Test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*0.9, linewidth=1, color='mediumaquamarine', linestyle='-.', label='test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.1, linewidth=2, color='mediumpurple', linestyle=':', label='test.')
    plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='mediumaquamarine', label='TALYS')
    plt.plot(energy_array, alice_spline(energy_array), linewidth=1, color='deepskyblue', linestyle='--', label='ALICE')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.2, linewidth=1, color='skyblue', label='Test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*0.9, linewidth=1, color='darkturquoise', linestyle='-.', label='test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.1, linewidth=2, color='cornflowerblue', linestyle=':', label='test.')


    plt.errorbar(energy_list[0:len(foil_list)], xs_list, xerr=[energy_min_unc_list[0:len(foil_list)], energy_plus_unc_list[0:len(foil_list)]], yerr=xs_unc_list, marker='D', markersize=5, linewidth=3, linestyle='', color='hotpink', label='This work')
    plt.xlabel('Beam energy (MeV)')
    plt.ylabel('Cross section (mb)')
    plt.title(reaction_product)
    plt.legend()
    plt.xlim(0,30)
    plt.show()




#_________________________________Running the code_______________________________________

# Need to change the function to R 

plot_xs('96NB', 41, 96, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05'])
plot_xs('90NB', 41, 90, ['Zr01', 'Zr02', 'Zr03', 'Zr04'])
# plot_xs('47SC', 21, 47, ['Ti01', 'Ti02', 'Ti03'])
# plot_xs('46SC', 21, 46, ['Ti01', 'Ti02', 'Ti03', 'Ti04'])
# plot_xs('48V', 23, 48, ['Ti01', 'Ti02', 'Ti03', 'Ti04'])
# plot_xs('52MN', 25, 52, ['Ni01', 'Ni02'])
# plot_xs('65NI', 28, 65, ['Ni01', 'Ni02', 'Ni03', 'Ni04'])
# plot_xs('57CO', 27, 57, ['Ni01', 'Ni02', 'Ni03'])
# plot_xs('60CU', 29, 60, ['Ni02', 'Ni03'])








