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
import curie as ci
import re
import requests
from foil_class_A0 import Foil_A0



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
    # A0_by_curie_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_curie.csv')
    # if os.path.exists(f'./Calculated_A0/{foil_name}_A0_by_hand.csv'):
    #     A0_by_hand_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_hand.csv')
    #     A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)
    # else:
    #     A0_concat_df = A0_by_curie_df
    A0_concat_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_hand.csv')


    A0_filtered_df = A0_concat_df[A0_concat_df['Isotope']==reaction_product]
    A0 = A0_filtered_df['A0'].iloc[0]
    A0_unc = A0_filtered_df['A0_unc'].iloc[0]

    foil = Foil_A0(foil_name, reaction_product, A0, A0_unc)
    foil.assign_areal_dens_w_unc_percent()
    foil.assign_molar_mass()
    foil.calculate_decay_constant_w_unc()
    foil.assign_beam_current_w_unc()
    foil.calculate_xs_w_unc_old()

    xs = foil.calc_xs
    xs_unc = foil.calc_xs_unc
    return xs, xs_unc #[mb]



def calc_xs_from_R(foil_name, reaction_product, reaction_product_w_state): #fungerer ikke med denne klassen
    #(foil_name, reaction_product, R, R_unc)
    R_by_curie_df = pd.read_csv(f'./Calculated_R/{foil_name}_R_by_curie.csv')
    print(R_by_curie_df)
    R_filtered_df = R_by_curie_df[R_by_curie_df['Isotope']==reaction_product_w_state]
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




def get_coh_data(path, filename):
    E_coh, xs_coh = np.loadtxt(path+filename, unpack = True)
    E_list = [0]+list(E_coh)
    xs_list = [0]+list(xs_coh)
    return E_list, xs_list


def get_empire_data(path, filename):
    E_empire, xs_empire = np.loadtxt(path+filename, unpack = True)
    E_list = [0]+list(E_empire)
    xs_list = [0]+list(xs_empire)
    return E_list, xs_list


def get_tendl_data(targetFoil, target, product, isomer):
    #example on how to run: get_tendl_data('Zr', 'Zr090', '039086')

    def tend_deuteron_url(targetFoil, target, product, fileEnding):
        return ('https://tendl.web.psi.ch/tendl_2023/deuteron_file/'
                + targetFoil + '/' + target
                + '/tables/residual/rp'
                + product + fileEnding)

    def tendl_file_ending(isomerLevel):
        return '.tot' if isomerLevel is None else '.L' + str(isomerLevel)

    def retrieve_tendl_data_from_url(url):
        tendl_data = requests.get(url).text.split("\n")[27:]  # skipping 27 first lines in tendl file
        tendl_data = np.genfromtxt(tendl_data)
        E = tendl_data[:, 0]
        Cs = tendl_data[:, 1]
        return E, Cs

    fileEnding = tendl_file_ending(isomer)
    url = tend_deuteron_url(targetFoil, target, product, fileEnding)
    E, xs = retrieve_tendl_data_from_url(url)
    return E, xs


def plot_xs(reaction_product, state, Z, A, foil_list, title, xs_type, exfor_manuel=False, write_csv=False, save_fig=False, use_R_to_calc_xs=True):

    plt.figure(figsize=(8, 6))

    target_material = foil_list[0][0:2]

    if target_material == 'Zr':
        target_material_big = 'ZR'
        target_A_list = [90, 91, 92, 94, 96]
    elif target_material == 'Ni':
        target_material_big = 'NI'
        target_A_list = [58, 60, 61, 62, 64]
    elif target_material == 'Ti':
        target_material_big = 'TI'
        target_A_list = [46, 47, 48, 49, 50]
    elif target_material == 'Fe':
        target_material_big = 'FE'
        target_A_list = [54, 56, 57, 58]


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
                   # 'Ti05': {'energy': 1.3297790055248617, 'min_unc': 0.9997790055248619, 'plus_unc': -0.3997790055248617},
                   'Ti05': {'energy': 1.3297790055248617, 'min_unc': 0.9997790055248619, 'plus_unc': 0},
                   'Fe01': {'energy': 48.2318781, 'min_unc': 0.8818780999999944, 'plus_unc': 0.9181219000000027},
                   'Zr06': {'energy': 47.63516853333334, 'min_unc': 0.8851685333333421, 'plus_unc': 0.914831466666655}, 
                   'Ti06': {'energy': 47.09795026666667, 'min_unc': 0.9479502666666662, 'plus_unc': 0.9520497333333253}, 
                   'Fe02': {'energy': 41.779047799999994, 'min_unc': 1.0290477999999936, 'plus_unc': 1.070952200000015}, 
                   'Zr07': {'energy': 41.110509566666664, 'min_unc': 1.060509566666667, 'plus_unc': 1.0394904333333415}, 
                   'Ti08': {'energy': 40.50592943333334, 'min_unc': 1.0559294333333398, 'plus_unc': 1.0440705666666545}, 
                   'Fe03': {'energy': 37.08828046666665, 'min_unc': 1.1382804666666502, 'plus_unc': 1.161719533333347}, 
                   'Zr08': {'energy': 36.351075, 'min_unc': 1.1010750000000016, 'plus_unc': 1.1989249999999956}, 
                   'Ti09': {'energy': 35.68120613333333, 'min_unc': 1.131206133333336, 'plus_unc': 1.1687938666666682}, 
                   'Fe04': {'energy': 31.872687066666668, 'min_unc': 1.2226870666666692, 'plus_unc': 1.277312933333338}, 
                   'Zr09': {'energy': 31.041554066666667, 'min_unc': 1.2915540666666665, 'plus_unc': 1.308445933333335}, 
                   'Ti10': {'energy': 30.28116373333334, 'min_unc': 1.3311637333333373, 'plus_unc': 1.3688362666666585}, 
                   'Fe05': {'energy': 25.879117466666663, 'min_unc': 1.4291174666666606, 'plus_unc': 1.470882533333338}, 
                   'Zr10': {'energy': 24.896101700000006, 'min_unc': 1.5461017000000048, 'plus_unc': 1.5538982999999966}, 
                   'Ti11': {'energy': 23.99016156666667, 'min_unc': 1.5401615666666686, 'plus_unc': 1.5598384333333293}
                   }

    energy_list = [energy_data[foil]['energy'] for foil in foil_list]
    energy_min_unc_list = [energy_data[foil]['min_unc'] for foil in foil_list]
    energy_plus_unc_list = [energy_data[foil]['plus_unc'] for foil in foil_list]

    xs_list = []
    xs_unc_list = []

    for foil in foil_list:
        if use_R_to_calc_xs==True:

            if state == 'ground_state':
                xs_state = 'g'
                reaction_product_w_state = reaction_product+xs_state
            elif state == 'isomeric_state':
                xs_state = 'm1'
                reaction_product_w_state = reaction_product+xs_state[-1]

            

            xs, xs_unc = calc_xs_from_R(foil, reaction_product, reaction_product_w_state)
            xs_list.append(xs)
            xs_unc_list.append(xs_unc)
        else:
            xs, xs_unc = calc_xs_from_A0(foil, reaction_product)
            xs_list.append(xs)
            xs_unc_list.append(xs_unc)


    #Write calculated xs and energy w unc to file
    if write_csv == True:
        csv_file_path = f'./Calculated_xs/{target_material}_dx_{reaction_product}_{state}_{xs_type}_xs.csv'
        with open(csv_file_path, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['Foil', 'Energy', 'Energy_min_unc', 'Energy_plus_unc', 'xs', 'xs_unc'])
            for i in range(len(foil_list)):
                csv_writer.writerow([f'{foil_list[i]}', f'{energy_list[i]}', f'{energy_min_unc_list[i]}', f'{energy_plus_unc_list[i]}', f'{xs_list[i]}', f'{xs_unc_list[i]}'])


    #Haleemas data
    E_86_haleema = [31.87, 36.34, 40.68, 47.23]
    dE_86_haleema = [1.5, 1.41, 1.14, 1.53]
    xs_86_haleema = [11.02, 29.80, 40.13, 32.29]
    d_xs_86_haleema = [0.84, 2.21, 2.98, 2.51]

    E_87_haleema = [18.43, 26.69, 25.10, 31.87, 36.34, 40.68, 47.23]
    dE_87_haleema = [1.49, 1.49, 1.41, 1.5, 1.41, 1.14, 1.53]
    xs_87_haleema = [10.08, 29.73, 26.60, 30.24, 24.19, 26.60, 36.55]
    d_xs_87_haleema = [1.64, 5.29, 1.97, 2.25, 1.79, 1.97, 2.83]

    E_88_haleema = [6.65, 12.20, 18.43, 26.69, 25.10, 31.87, 36.34, 40.68, 47.23]
    dE_88_haleema = [1.78, 1.58, 1.49, 1.49, 1.41, 1.5, 1.41, 1.14, 1.53]
    xs_88_haleema = [0.4, 6.25, 8.18, 12.57, 11.67, 15.43, 25.10, 40.77, 58.23]
    d_xs_88_haleema = [0.07, 0.89, 1.27, 2.26, 0.88, 1.18, 1.90, 3.06, 4.57]


    #Get data from talys
    talys_file = f'rp0{Z}0{A}.tot'
    E_talys, xs_talys = get_talys_data(talys_file, foil_list[0])
    talys_spline = splrep(E_talys, xs_talys)
    energy_array = np.linspace(0,50,500)


    #Get data from alice
    E_alice, xs_alice = get_alice_data(target_material, Z, A)
    E_alice = [0]+E_alice
    xs_alice = [0]+xs_alice
    alice_spline = PchipInterpolator(E_alice, xs_alice)


    #Get data from coh
    E_coh_list = []
    xs_coh_list = []
    abundances_coh_list = []

    el = ci.Element(target_material)
    abundances_df = el.abundances
   
    if state == 'ground_state':
        coh_state = 'G'
    elif state == 'isomeric_state':
        coh_state = 'M'

    for i in range(len(target_A_list)):
        coh_path = f'./coh_calculations/{target_A_list[i]}{target_material}/'
        if state == 'isomeric_state':
            product = reaction_product[:-1]
        else:
            product = reaction_product

        coh_file = f'0{Z}-0{A}{product[2:]}{coh_state}_coh.txt'
        # Check if the file exists:
        if os.path.exists(coh_path + coh_file):
            E_coh_partial, xs_coh_partial = get_coh_data(coh_path, coh_file)
            E_coh_list.append(E_coh_partial)
            xs_coh_list.append(xs_coh_partial)

            isotope = f'{target_A_list[i]}{target_material_big}'
            ab = abundances_df[abundances_df['isotope'] == isotope]['abundance'].values[0]
            abundances_coh_list.append(ab)

    xs_coh = [0] * len(xs_coh_list[0])
    # Loop through each sublist and multiply it with the corresponding abundance, then sum them up
    for i, sublist in enumerate(xs_coh_list):
        abundance = abundances_coh_list[i]
        weighted_sublist = [value * abundance/100 for value in sublist]
        xs_coh = [x + y for x, y in zip(xs_coh, weighted_sublist)]

    E_coh = E_coh_list[0]

    coh_spline = PchipInterpolator(E_coh, xs_coh)


    

    #Get data from empire
    E_empire_list = []
    xs_empire_list = []
    abundances_empire_list = []

    el = ci.Element(target_material)
    abundances_df = el.abundances
   
    if state == 'ground_state':
        empire_state = 'G'
    elif state == 'isomeric_state':
        empire_state = 'M'

    for i in range(len(target_A_list)):
        empire_path = f'./empire_calculations/{target_material}/{target_A_list[i]}{target_material}/xs_plots/'
        if state == 'isomeric_state':
            product = reaction_product[:-1]
        else:
            product = reaction_product

        empire_file = f'{Z}-{product[2:]}-{A}{empire_state}_empire.txt'
        
        if os.path.exists(empire_path + empire_file):
            E_empire_partial, xs_empire_partial = get_empire_data(empire_path, empire_file)

            E_empire_list.append(E_empire_partial)
            xs_empire_list.append(xs_empire_partial)

            isotope = f'{target_A_list[i]}{target_material_big}'
            ab = abundances_df[abundances_df['isotope'] == isotope]['abundance'].values[0]
            abundances_empire_list.append(ab)

    xs_empire = [0] * len(xs_empire_list[0])
    # Loop through each sublist and multiply it with the corresponding abundance, then sum them up
    for i, sublist in enumerate(xs_empire_list):
        abundance = abundances_empire_list[i]
        weighted_sublist = [value * abundance/100 for value in sublist]
        xs_empire = [x + y for x, y in zip(xs_empire, weighted_sublist)]

    E_empire = E_empire_list[0]
    empire_spline = PchipInterpolator(E_empire, xs_empire)



    #Get data from tendl
    E_tendl_list = []
    xs_tendl_list = []
    abundances_tendl_list = []

    el = ci.Element(target_material)
    abundances_df = el.abundances

    if state == 'ground_state':
        isomerLevel = None
    elif state == 'isomeric_state':
        isomerLevel = '00'

    for i in range(len(target_A_list)):   
        try:
            E_tendl_partial, xs_tendl_partial = get_tendl_data(target_material, f'{target_material}0{target_A_list[i]}', f'0{Z}0{A}', isomerLevel)
            
            if E_tendl_partial is not None and xs_tendl_partial is not None:
                E_tendl_list.append(E_tendl_partial)
                xs_tendl_list.append(xs_tendl_partial)

                isotope = f'{target_A_list[i]}{target_material_big}'
                ab = abundances_df[abundances_df['isotope'] == isotope]['abundance'].values[0]
                abundances_tendl_list.append(ab)
            else:
                print(f'Skipping data retrieval for A={target_A_list[i]} due to an issue with the URL or data format.')
        except Exception as e:
            print(f'Unable to retrieve tendl data for A={target_A_list[i]}\n{str(e)}')
            continue

    xs_tendl = [0] * len(xs_tendl_list[0])

    # Loop through each sublist and multiply it with the corresponding abundance, then sum them up
    for i, sublist in enumerate(xs_tendl_list):
        abundance = abundances_tendl_list[i]
        weighted_sublist = [value * abundance/100 for value in sublist]
        xs_tendl = [x + y for x, y in zip(xs_tendl, weighted_sublist)]

    E_tendl = E_tendl_list[0]
    tendl_spline = PchipInterpolator(E_tendl, xs_tendl)



    #Get data from exfor
    markers = ['.', '*', 'v', '^', '+', '<', '>', 's', 'h',     '.', '*', 'v', '^', '+', '<', '>', 's', 'h', '.', '*', 'v', '^', '+', '<', '>', 's', 'h',     '.', '*', 'v', '^', '+', '<', '>', 's', 'h']
    grey_colors = ['dimgrey', 'darkgrey', 'lightgrey', 'silver', 'k', 'dimgrey', 'darkgrey', 'lightgrey', 'silver',     'silver', 'lightgrey', 'darkgrey', 'dimgrey', 'k','silver', 'lightgrey', 'dimgrey', 'darkgrey', 'darkgrey', 'lightgrey', 'silver', 'k', 'dimgrey', 'darkgrey', 'lightgrey', 'silver',     'silver', 'lightgrey', 'darkgrey', 'dimgrey', 'k','silver', 'lightgrey', 'dimgrey', 'darkgrey', 'dimgrey']
    # ref_numbers = ['[66]', '[57]', '[62]', '[58]', '[59]', '[64]']
    ref_numbers = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']

    if exfor_manuel == False:
        product = reaction_product[2:]+'-'+reaction_product[0:2]

        if foil_list[0][0:2]=='Ti':
            target = 'TI-0'
        elif foil_list[0][0:2]=='Zr':
            target = 'ZR-0'
        elif foil_list[0][0:2]=='Ni':
            target = 'NI-0'
        elif foil_list[0][0:2]=='Fe':
            target = 'FE-0'

        try:
            exfor_dict = get_exfor_data(target, 'D,*', product)
             #Plot results
            
            k=0
            for index in exfor_dict:
                # Need to set up list of marker sizes to iterate over with k
                # print(k)
                if exfor_dict[index][2].shape[1] == 4:
                    plt.errorbar(exfor_dict[index][2][:,0],1E3*exfor_dict[index][2][:,1], xerr=exfor_dict[index][2][:,2], yerr=1E3*exfor_dict[index][2][:,3], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')'+f' {ref_numbers[k]}', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])

                elif exfor_dict[index][2].shape[1] == 3:
                    plt.errorbar(exfor_dict[index][2][:,0],1E3*exfor_dict[index][2][:,1], yerr=1E3*exfor_dict[index][2][:,2], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')'+f' {ref_numbers[k]}', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])

                elif exfor_dict[index][2].shape[1] >= 5:
                    print('WARNING: Plotting',str(exfor_dict[index][2].shape[1])+'-column EXFOR data retrieved for subentry', exfor_dict[index][3]+', please make sure data look reasonable - column formatting is inconsistent for >4 columns.')
                    plt.errorbar(exfor_dict[index][2][:,4],exfor_dict[index][2][:,5], xerr=exfor_dict[index][2][:,0], yerr=exfor_dict[index][2][:,6], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')'+f' {ref_numbers[k]}', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])
                k=k+1

        except UnboundLocalError:
            print('Warning: Not able to plot exfor data!')


    elif exfor_manuel == True:
        path = f'./exfor/{foil_list[0][0:2]}/{reaction_product}_{xs_type}'
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
            label = f"{author_name} ({year}) {ref_numbers[i]}"


            plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)

    cum_xs_87Yg = np.array([28.866671386714938, 10.00821348815171, 39.37983486276633, 26.90301244292157, 24.350239464668082, 31.77747524552684, 29.427202521178952])
    ind_xs_87Yg = np.array([6.96835117734811, 3.1001232016703795, 7.343434769274522, 6.131376778114652, 5.958852822320928, 7.061281585995764, 6.77368726701053])
    # cum_xs_87Ym = np.array([21.694067175519578, 6.700992387276112, 18.138671979583105, 16.152673100433443, 15.527035461932755, 20.14399331442515, 16.04748704075385])
    cum_xs_87Ym = np.array([21.694067175519578, 6.700992387276112, 32.503102, 21.123407, 18.701477, 25.021261, 23.051401])
    energy = np.array([26.395312639999993, 19.241718260000003, 47.63516853333334, 41.110509566666664, 36.351075, 31.041554066666667, 24.896101700000006])




    diff_xs = cum_xs_87Yg-cum_xs_87Ym


    # plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='gold', label='TALYS')
    # plt.plot(energy_array, splev(energy_array, alice_spline), linewidth=1, color='plum', linestyle='--', label='ALICE')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.2, linewidth=1, color='skyblue', label='Test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*0.9, linewidth=1, color='mediumaquamarine', linestyle='-.', label='test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.1, linewidth=2, color='mediumpurple', linestyle=':', label='test.')

    if reaction_product == '86Y':
        plt.errorbar(E_86_haleema, xs_86_haleema, xerr=dE_86_haleema, yerr=d_xs_86_haleema, linewidth = 1, capsize = 1, marker='s', linestyle='None', color='black', label= 'Zaneb (2018) [24]')
    if reaction_product == '87Y':
        plt.errorbar(E_87_haleema, xs_87_haleema, xerr=dE_87_haleema, yerr=d_xs_87_haleema, linewidth = 1, capsize = 1, marker='s', linestyle='None', color='black', label= 'Zaneb (2018) [24]')
    if reaction_product == '88Y':
        plt.errorbar(E_88_haleema, xs_88_haleema, xerr=dE_88_haleema, yerr=d_xs_88_haleema, linewidth = 1, capsize = 1, marker='s', linestyle='None', color='black', label= 'Zaneb (2018) [24]')

    plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='mediumseagreen', label='TALYS-2.0')
    plt.plot(energy_array, alice_spline(energy_array), linewidth=1, color='mediumblue', linestyle=(0,(5,5)), label='ALICE-2020')
    plt.plot(energy_array, coh_spline(energy_array), linewidth=1, color='olivedrab', linestyle=(0,(12,10)), label='CoH-3.5.3')
    plt.plot(energy_array, empire_spline(energy_array), linewidth=1.5, color='deepskyblue', linestyle=(0,(5,4,1,4,1,4)), label='EMPIRE 3.2.3')
    plt.plot(energy_array, tendl_spline(energy_array), linewidth=2, color='cornflowerblue', linestyle=':', label='TENDL-2023')

    # plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='mediumaquamarine', label='TALYS-2.0')
    # plt.plot(energy_array, alice_spline(energy_array), linewidth=1, color='deepskyblue', linestyle=(0,(5,1)), label='ALICE-2020')
    # plt.plot(energy_array, coh_spline(energy_array), linewidth=1, color='skyblue', linestyle=(5,(10,4)), label='CoH-3.5.3')
    # plt.plot(energy_array, empire_spline(energy_array), linewidth=1, color='darkturquoise', linestyle=(0,(5,2,1,2,1,2)), label='EMPIRE 3.2.3')
    # plt.plot(energy_array, tendl_spline(energy_array), linewidth=2, color='cornflowerblue', linestyle=':', label='TENDL-2023')



    # plt.errorbar(energy_list[0:len(foil_list)], xs_list, xerr=[energy_min_unc_list[0:len(foil_list)], energy_plus_unc_list[0:len(foil_list)]], yerr=xs_unc_list, marker='D', markersize=5, linewidth=3, linestyle='', color='hotpink', label='This work')
    plt.errorbar(energy_list, xs_list, xerr=[energy_min_unc_list, energy_plus_unc_list], yerr=xs_unc_list, marker='D', markersize=4, linewidth=2, linestyle='', color='hotpink', label='This work')

    plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
    plt.ylabel('Cross Section (mb)', fontsize=14)
    plt.title(f'{title}{xs_type}', fontsize=14)
    # plt.legend(fontsize=12, ncol=2)
    plt.legend(fontsize=12)
    # plt.legend(fontsize=4)
    plt.ylim((0))
    plt.xlim(0,50)
    if save_fig==True:
        plt.savefig(f'./Figures/PhD/xs_plots/{target_material}/{target_material}_dx_{reaction_product}_{state}_{xs_type}.pdf', dpi=600)

    # plt.plot(energy, diff_xs, 'bo')
    plt.show()
    









#_________________________________Running the code______________________________________________________________________________________________________________

# #_______________________________Zr-foils________________________________________________
# plot_xs('88Y', 'ground_state', 39, 88, ['Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{88}$Y - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('88ZR', 'ground_state', 40, 88, ['Zr06', 'Zr07'], r'$^{nat}$Zr(d,x)$^{88}$Zr - ', 'Independent', save_fig=True, write_csv=False)
# plot_xs('88NB', 'ground_state', 41, 88, ['Zr06', 'Zr07'], r'$^{nat}$Zr(d,x)$^{88}$Nb - ', 'Independent', save_fig=True, write_csv=False)

# plot_xs('90NB', 'ground_state', 41, 90, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{90}$Nb - ', 'Independent', save_fig=True, write_csv=False)

# plot_xs('92NBm', 'isomeric_state', 41, 92, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{92m}$Nb - ', 'Independent', exfor_manuel=True, save_fig=True, write_csv=False)

# plot_xs('95NBm', 'isomeric_state', 41, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{95m}$Nb - ', 'Independent', save_fig=True, exfor_manuel=True, write_csv=False)
# plot_xs('95NB', 'ground_state', 41, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{95}$Nb - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)

# plot_xs('96NB', 'ground_state', 41, 96, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{96}$Nb - ', 'Independent', save_fig=True, write_csv=False)





# # _______________________________Ni-foils________________________________________________
# plot_xs('54MN', 'ground_state', 25, 54, ['Ni01'], r'$^{nat}$Ni(d,x)$^{54}$Mn - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)

# plot_xs('57CO', 'ground_state', 27, 57, ['Ni01', 'Ni02', 'Ni03'], r'$^{nat}$Ni(d,x)$^{57}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True) 

# plot_xs('65NI', 'ground_state', 28, 65, ['Ni01', 'Ni02', 'Ni03', 'Ni04'], r'$^{nat}$Ni(d,x)$^{65}$Ni - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)





# # _______________________________Fe-foils________________________________________________
# plot_xs('54MN', 'ground_state', 25, 54, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{54}$Mn - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('55CO', 'ground_state', 27, 55, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{55}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('57CO', 'ground_state', 27, 57, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{57}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('58CO', 'ground_state', 27, 58, ['Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{58}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
































































# #_______________________________Zr-foils________________________________________________
# # plot_xs('86Y', 'ground_state', 39, 86, ['Zr06', 'Zr07', 'Zr08', 'Zr09'], r'$^{nat}$Zr(d,x)$^{86}$Y - ', 'Independent', save_fig=True, write_csv=False)
# plot_xs('86Ym', 'isomeric_state', 39, 86, ['Zr06', 'Zr07', 'Zr08', 'Zr09'], r'$^{nat}$Zr(d,x)$^{86m}$Y - ', 'Independent', save_fig=False, write_csv=False, exfor_manuel=True)

# plot_xs('87Y', 'ground_state', 39, 87, ['Zr01', 'Zr02', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{87}$Y - ', 'Independent', save_fig=False, write_csv=True, exfor_manuel=True)
# # # plot_xs('87Ym', 'isomeric_state', 39, 87, ['Zr01', 'Zr02', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{87m}$Y - ', 'Cumulative', save_fig=True, write_csv=True, exfor_manuel=False)

# # # plot_xs('87Ym', 'isomeric_state', 39, 87, ['Zr01', 'Zr02', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{87m}$Y - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=False)
# # # plot_xs('87ZR', 'ground_state', 40, 87, ['Zr01', 'Zr02', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{87}$Zr - ', 'Cumulative', save_fig=True, write_csv=False)


# # plot_xs('88Y', 'ground_state', 39, 88, ['Zr01', 'Zr02', 'Zr03', 'Zr04'], r'$^{nat}$Zr(d,x)$^{88}$Y - ', 'Cumulative', save_fig=True, write_csv=False)
# plot_xs('88Y', 'ground_state', 39, 88, ['Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{88}$Y - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('88ZR', 'ground_state', 40, 88, ['Zr06', 'Zr07'], r'$^{nat}$Zr(d,x)$^{88}$Zr - ', 'Independent', save_fig=True, write_csv=False)
# # plot_xs('88ZR', 'ground_state', 40, 88, ['Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{88}$Zr - ', 'Cumulative', save_fig=True, write_csv=False)
# plot_xs('88NB', 'ground_state', 41, 88, ['Zr06', 'Zr07'], r'$^{nat}$Zr(d,x)$^{88}$Nb - ', 'Independent', save_fig=True, write_csv=False)

# # # plot_xs('88Y', 'ground_state', 39, 88, ['Zr01', 'Zr02', 'Zr03', 'Zr04'], r'$^{nat}$Zr(d,x)$^{88}$Y - ', 'Cumulative', save_fig=True, write_csv=False)
# # # plot_xs('88Y', 'ground_state', 39, 88, ['Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{88}$Y - ', 'Independent', save_fig=True, write_csv=False)
# # # plot_xs('88ZR', 'ground_state', 40, 88, ['Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{88}$Zr - ', 'Cumulative', save_fig=True, write_csv=False)

# # # plot_xs('89NBm', 'isomeric_state', 41, 89, ['Zr01', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{89m}$Nb - ', 'Independent', save_fig=True, write_csv=False)
# # # plot_xs('89NB', 'ground_state', 41, 89, ['Zr01', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{89}$Nb - ', 'Cumulative', save_fig=True, write_csv=False)
# # # plot_xs('89ZR', 'ground_state', 40, 89, ['Zr01', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{89}$Zr - ', 'Independent', save_fig=True, write_csv=False)
# # plot_xs('89ZR', 'ground_state', 40, 89, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{89}$Zr - ', 'Cumulative', save_fig=True, write_csv=False)


# plot_xs('90NB', 'ground_state', 41, 90, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{90}$Nb - ', 'Independent', save_fig=True, write_csv=False)

# # plot_xs('90Ym', 'isomeric_state', 39, 90, ['Zr01', 'Zr03', 'Zr06', 'Zr07', 'Zr08', 'Zr10'], r'$^{nat}$Zr(d,x)$^{90m}$Y - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=True)

# plot_xs('92NBm', 'isomeric_state', 41, 92, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{92m}$Nb - ', 'Independent', exfor_manuel=True, save_fig=True, write_csv=False)

# # plot_xs('95ZR', 'ground_state', 40, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{95}$Zr - ', 'Cumulative', save_fig=True, write_csv=False)
# plot_xs('95NBm', 'isomeric_state', 41, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{95m}$Nb - ', 'Independent', save_fig=True, exfor_manuel=True, write_csv=False)
# plot_xs('95NB', 'ground_state', 41, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{95}$Nb - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)

# plot_xs('96NB', 'ground_state', 41, 96, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{96}$Nb - ', 'Independent', save_fig=False, write_csv=False)





# # _______________________________Ni-foils________________________________________________
# plot_xs('52MN', 'ground_state', 25, 52, ['Ni01', 'Ni02'], r'$^{nat}$Ni(d,x)$^{52}$Mn - ', 'Cumulative', save_fig=False, write_csv=False, exfor_manuel=True)
# plot_xs('54MN', 'ground_state', 25, 54, ['Ni01'], r'$^{nat}$Ni(d,x)$^{54}$Mn - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('55CO', 'ground_state', 27, 55, ['Ni01', 'Ni02', 'Ni03'], r'$^{nat}$Ni(d,x)$^{55}$Co - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('56CO', 'ground_state', 27, 56, ['Ni01', 'Ni02', 'Ni03', 'Ni04'], r'$^{nat}$Ni(d,x)$^{56}$Co - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('57NI', 'ground_state', 28, 57, ['Ni01', 'Ni02', 'Ni03'], r'$^{nat}$Ni(d,x)$^{57}$Ni - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('57CO', 'ground_state', 27, 57, ['Ni01', 'Ni02', 'Ni03'], r'$^{nat}$Ni(d,x)$^{57}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True) #Har resultater for Ni05, men stor unc
# plot_xs('57CO', 'ground_state', 27, 57, ['Ni04'], r'$^{nat}$Ni(d,x)$^{57}$Co - ', 'Cumulative', save_fig=False, write_csv=False, exfor_manuel=False) #Har resultater for Ni05, men stor unc

# plot_xs('58CO', 'ground_state', 27, 58, ['Ni01', 'Ni02', 'Ni03', 'Ni04'], r'$^{nat}$Ni(d,x)$^{58}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True) #Har resultater for Ni05, men stor unc
# plot_xs('60CO', 'ground_state', 27, 60, ['Ni01',  'Ni03', 'Ni04'], r'$^{nat}$Ni(d,x)$^{60}$Co - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=True) 
# plot_xs('61CU', 'ground_state', 29, 61, ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05'], r'$^{nat}$Ni(d,x)$^{61}$Cu - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('64CU', 'ground_state', 29, 64, ['Ni02', 'Ni03', 'Ni04'], r'$^{nat}$Ni(d,x)$^{64}$Cu - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('65NI', 'ground_state', 28, 65, ['Ni01', 'Ni02', 'Ni03', 'Ni04'], r'$^{nat}$Ni(d,x)$^{65}$Ni - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)





# # _______________________________Ti-foils________________________________________________

# plot_xs('44SC', 'ground_state', 21, 44, ['Ti06'], r'$^{nat}$Ti(d,x)$^{44}$Sc - ', 'Independent', save_fig=False, exfor_manuel=True)
# plot_xs('44SCm', 'isomeric_state', 21, 44, ['Ti06'], r'$^{nat}$Ti(d,x)$^{44m}$Sc - ', 'Cumulative', save_fig=False, exfor_manuel=True)
# plot_xs('46SC', 'ground_state', 21, 46, ['Ti01', 'Ti02', 'Ti03', 'Ti04', 'Ti05', 'Ti06', 'Ti08', 'Ti09', 'Ti10', 'Ti11'], r'$^{nat}$Ti(d,x)$^{46}$Sc - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('47SC', 'ground_state', 21, 47, ['Ti01', 'Ti02', 'Ti03', 'Ti06', 'Ti08', 'Ti09', 'Ti10', 'Ti11'], r'$^{nat}$Ti(d,x)$^{47}$Sc - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('48SC', 'ground_state', 21, 48, ['Ti06', 'Ti08', 'Ti09', 'Ti10'], r'$^{nat}$Ti(d,x)$^{48}$Sc - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('48V', 'ground_state', 23, 48, ['Ti01', 'Ti02', 'Ti03', 'Ti04', 'Ti05', 'Ti06', 'Ti08', 'Ti09', 'Ti10', 'Ti11'], r'$^{nat}$Ti(d,x)$^{48}$V - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)



# # _______________________________Fe-foils________________________________________________

# # plot_xs('48V', 'ground_state', 23, 48, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{48}$V - ', 'Cumulative', save_fig=False, write_csv=False, exfor_manuel=True)
# # plot_xs('52MN', 'ground_state', 25, 52, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{52}$Mn - ', 'Cumulative', save_fig=False, write_csv=False, exfor_manuel=True)
# plot_xs('54MN', 'ground_state', 25, 54, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{54}$Mn - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('55CO', 'ground_state', 27, 55, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{55}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# # plot_xs('56CO', 'ground_state', 27, 56, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{56}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('57CO', 'ground_state', 27, 57, ['Fe01', 'Fe02', 'Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{57}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)
# plot_xs('58CO', 'ground_state', 27, 58, ['Fe03', 'Fe04', 'Fe05'], r'$^{nat}$Fe(d,x)$^{58}$Co - ', 'Independent', save_fig=True, write_csv=False, exfor_manuel=True)















