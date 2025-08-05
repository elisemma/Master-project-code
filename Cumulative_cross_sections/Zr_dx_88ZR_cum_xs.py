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


df_88Y_cum = pd.read_csv('../Calculated_xs/Zr_dx_88Y_ground_state_Cumulative_xs.csv')
df_88Y_ind = pd.read_csv('../Calculated_xs/Zr_dx_88Y_ground_state_Independent_xs.csv')
df_88Zr_cum = pd.read_csv('../Calculated_xs/Zr_dx_88ZR_ground_state_Cumulative_xs.csv')
df_88Zr_ind = pd.read_csv('../Calculated_xs/Zr_dx_88ZR_ground_state_Independent_xs.csv')
df_88Nb_ind = pd.read_csv('../Calculated_xs/Zr_dx_88NB_ground_state_Independent_xs.csv')


# 88-chain: I want the cumualtive xs for 88Y. I have the cumulative xs for 88Y in foil 1,2,3,4
# I need to add with branching ratios for foil 6,7,8,9,10
# For foil 8,9,10: I know the cumulative xs for 88Zr and can just add with branchingratio
# For foil 6,7: I first need to find the cumuative xs for 88Zr by using the xs for 88Nb

# 88Nb decays to 88Zr (eps = 100%)
# 88Zr decays to 88Y (eps = 100%)

#Since the branching ratios are 100%, I can just add the cross sections:


df_88Nb_ind = df_88Nb_ind[::-1]
df_88Zr_ind = df_88Zr_ind[::-1]

# Cumulative xs for 88Zr:
print('df_88Zr_cum before adding anything:')
print(df_88Zr_cum)

for row1, row2 in zip(df_88Nb_ind.itertuples(), df_88Zr_ind.itertuples()):

    combined_xs = row1.xs + row2.xs 
    combined_xs_unc = np.sqrt(row1.xs_unc**2 + row2.xs_unc**2)
    print(row1.Foil, row2.Foil)
    new_row = pd.DataFrame({
        'Foil': row1.Foil,
        'Energy': row1.Energy,
        'Energy_min_unc': row1.Energy_min_unc,
        'Energy_plus_unc': row1.Energy_plus_unc,
        'xs': combined_xs,
        'xs_unc': combined_xs_unc
    }, index=[0])

    # df_88Zr_cum = df_88Zr_cum._append(new_row, ignore_index=True)

    df_88Zr_cum = pd.concat([new_row, df_88Zr_cum], ignore_index=True)

print('\ndf_88Zr_cum after adding 88Nb_ind of 88_Zr_ind:')
print(df_88Zr_cum)



print('\n')


# Cumulative xs for 88Y:
print('\ndf_88Y_cum before adding anything:')
print(df_88Y_cum)

for row1, row2 in zip(df_88Zr_cum.itertuples(), df_88Y_ind.itertuples()):

    combined_xs = row1.xs + row2.xs 
    combined_xs_unc = np.sqrt(row1.xs_unc**2 + row2.xs_unc**2)
    print(row1.Foil, row2.Foil)
    new_row = {
        'Foil': row1.Foil,
        'Energy': row1.Energy,
        'Energy_min_unc': row1.Energy_min_unc,
        'Energy_plus_unc': row1.Energy_plus_unc,
        'xs': combined_xs,
        'xs_unc': combined_xs_unc
    }

    df_88Y_cum = df_88Y_cum._append(new_row, ignore_index=True)


print('\ndf_88Y_cum after adding 88Y_ind of 88_Zr_cum:')
print(df_88Y_cum)







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
    A0_by_curie_df = pd.read_csv(f'../Calculated_A0/{foil_name}_A0_by_curie.csv')
    if os.path.exists(f'../Calculated_A0/{foil_name}_A0_by_hand.csv'):
        A0_by_hand_df = pd.read_csv(f'../Calculated_A0/{foil_name}_A0_by_hand.csv')
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



def calc_xs_from_R(foil_name, reaction_product, reaction_product_w_state): #fungerer ikke med denne klassen
    #(foil_name, reaction_product, R, R_unc)
    R_by_curie_df = pd.read_csv(f'../Calculated_R/{foil_name}_R_by_curie.csv')
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
    E_talys, xs_talys= np.loadtxt(f'../talys_calculations/talys_{foil_name[0:2]}/talys/'+filename, unpack = True, skiprows = 24)
    E_feeding, xs_feeding= np.loadtxt(f'../talys_calculations/talys_{foil_name[0:2]}/talys/rp041088.tot', unpack = True, skiprows = 24)

    E_list = [0]+list(E_talys)
    xs_ind = [0]+list(xs_talys)
    xs_feeding = [0]+list(xs_feeding)

    xs_list = [x+y for x,y in zip(xs_ind, xs_feeding)]
    return E_list, xs_list



def get_alice_data(target_material, Z, A):
    # Initialize lists to store energy and cross-section data
    file_path = f'../alice_calculations/plot_nat{target_material}_dx'
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
                        cross_sections.append(float(data[3])) #[3] er total xs, [5] er gs, [7] er isomer 1
                except ValueError:
                    pass  # Skip lines that can't be converted to integers
    
    return energies, cross_sections



def get_empire_data(path, filename):
    E_empire, xs_empire = np.loadtxt(path+filename, unpack = True)
    E_list = [0]+list(E_empire)
    xs_ind = [0]+list(xs_empire)

    if os.path.exists(path + '41-Nb-88G_empire.txt'):
        E_feeding, xs_feeding = np.loadtxt(path + '41-Nb-88G_empire.txt', unpack = True)
        xs_feeding = [0]+list(xs_feeding)
        xs_list = [x+y for x,y in zip(xs_ind, xs_feeding)]
    else:
        xs_list = xs_ind

    return E_list, xs_list



def get_coh_data(path, filename):
    E_coh, xs_coh = np.loadtxt(path+filename, unpack = True)
    E_list = [0]+list(E_coh)
    xs_ind = [0]+list(xs_coh)
    xs_list = xs_ind

    if os.path.exists(path + '041-088NbG_coh.txt'):
        E_feeding, xs_feeding = np.loadtxt(path+'041-088NbG_coh.txt', unpack = True)
        xs_feeding = [0]+list(xs_feeding)
        xs_list = [x+y for x,y in zip(xs_list, xs_feeding)]
    if os.path.exists(path + '040-088ZrM_coh.txt'):
        E_feeding2, xs_feeding2 = np.loadtxt(path+'040-088ZrM_coh.txt', unpack = True)
        xs_feeding2 = [0]+list(xs_feeding2)
        xs_list = [x+y for x,y in zip(xs_list, xs_feeding2)]

    return E_list, xs_list


def get_tendl_data(targetFoil, target, product, isomerLevel):
    #example on how to run: get_tendl_data('Zr', 'Zr090', '039086', None)

    def tend_deuteron_url(targetFoil, target, product, fileEnding):
        return ('https://tendl.web.psi.ch/tendl_2023/deuteron_file/'
                + targetFoil + '/' + target
                + '/tables/residual/rp'
                + product + fileEnding)

    def tendl_file_ending(isomerLevel):
        return '.tot' if isomerLevel is None else '.L' + str(isomerLevel)

    def retrieve_tendl_data_from_url(url):
        tendl_data = requests.get(url).text.split("\n")[27:]  #skipping 27 first lines in tendl file
        tendl_data = np.genfromtxt(tendl_data)
        E = tendl_data[:, 0]
        Cs = tendl_data[:, 1]
        return E, Cs

    fileEnding = tendl_file_ending(isomerLevel)
    url = tend_deuteron_url(targetFoil, target, product, fileEnding)
    E, xs = retrieve_tendl_data_from_url(url)
    xs_tot = xs


    try:
        fileEnding_feeding1 = tendl_file_ending(None)
        url_feeding1 = tend_deuteron_url(targetFoil, target, '041088', fileEnding_feeding1)
        E_feeding1, xs_feeding1 = retrieve_tendl_data_from_url(url_feeding1)
        xs_tot = [x+y for x,y in zip(xs, xs_feeding1)]
    except:
        print(f'No feeding for {target}')

    return E, xs_tot


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


    #Write calculated xs and energy w unc to file
    if write_csv == True:
        csv_file_path = f'../Calculated_xs/{target_material}_dx_{reaction_product}_{state}_{xs_type}_xs.csv'
        with open(csv_file_path, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['Foil', 'Energy', 'Energy_min_unc', 'Energy_plus_unc', 'xs', 'xs_unc'])
            # for i in range(len(foil_list)):
            #     csv_writer.writerow([f'{foil_list[i]}', f'{energy_list[i]}', f'{energy_min_unc_list[i]}', f'{energy_plus_unc_list[i]}', f'{xs_list[i]}', f'{xs_unc_list[i]}'])
            for row in df_88Zr_cum.itertuples():
                csv_writer.writerow([f'{row.Foil}', f'{row.Energy}', f'{row.Energy_min_unc}', f'{row.Energy_plus_unc}', f'{row.xs}', f'{row.xs_unc}'])



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
    E_alice_feed, xs_alice_feed = get_alice_data(target_material, 41, 88)
    xs_alice_list = [0]+xs_alice[0:4]+[x+y for x,y in zip(xs_alice[4:],xs_alice_feed)]
    alice_spline = PchipInterpolator(E_alice, xs_alice_list)
    
    
 
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
        empire_path = f'../empire_calculations/{target_material}/{target_A_list[i]}{target_material}/xs_plots/'
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
        coh_path = f'../coh_calculations/{target_A_list[i]}{target_material}/'
        if state == 'isomeric_state':
            product = reaction_product[:-1]
        else:
            product = reaction_product

        coh_file = f'0{Z}-0{A}{product[2:]}{coh_state}_coh.txt'

        # coh_file_isom = f'0{Z}-0{A}{product[2:]}M_coh.txt'

        # Check if the file exists:
        if os.path.exists(coh_path + coh_file):
            E_coh_partial, xs_coh_partial = get_coh_data(coh_path, coh_file)

            # E_coh_partial_isom, xs_coh_partial_isom = get_coh_data(coh_path, coh_file_isom)
            # for j in range(len(E_coh_partial)):
            #     xs_coh_partial[j]+=xs_coh_partial_isom[j]

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
    markers = ['.', '*', 'v', '^', '+', '<', '>', 's', 'h',     '.', '*', 'v', '^', '+', '<', '>', 's', 'h']
    grey_colors = ['dimgrey', 'darkgrey', 'lightgrey', 'silver', 'k', 'dimgrey', 'darkgrey', 'lightgrey', 'silver',     'silver', 'lightgrey', 'darkgrey', 'dimgrey', 'k','silver', 'lightgrey', 'dimgrey', 'darkgrey']

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
                print(exfor_dict[index][0])
                if exfor_dict[index][2].shape[1] == 4:
                    if exfor_dict[index][0] == 'Vysotskij':
                        plt.errorbar(exfor_dict[index][2][:,0],exfor_dict[index][2][:,1], xerr=exfor_dict[index][2][:,2], yerr=exfor_dict[index][2][:,3], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])
                    else:
                        plt.errorbar(exfor_dict[index][2][:,0],1E3*exfor_dict[index][2][:,1], xerr=exfor_dict[index][2][:,2], yerr=1E3*exfor_dict[index][2][:,3], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])

                elif exfor_dict[index][2].shape[1] == 3:
                    if exfor_dict[index][0] == 'Vysotskij':
                        plt.errorbar(exfor_dict[index][2][:,0],exfor_dict[index][2][:,1], yerr=exfor_dict[index][2][:,2], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])
                    else:
                        plt.errorbar(exfor_dict[index][2][:,0],1E3*exfor_dict[index][2][:,1], yerr=1E3*exfor_dict[index][2][:,2], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])
                elif exfor_dict[index][2].shape[1] >= 5:
                    if exfor_dict[index][0] == 'Vysotskij':

                        print('WARNING: Plotting',str(exfor_dict[index][2].shape[1])+'-column EXFOR data retrieved for subentry', exfor_dict[index][3]+', please make sure data look reasonable - column formatting is inconsistent for >4 columns.')
                        plt.errorbar(exfor_dict[index][2][:,4],1E-3*exfor_dict[index][2][:,5], xerr=exfor_dict[index][2][:,0], yerr=1E-3*exfor_dict[index][2][:,6], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])
                    else:
                        print('WARNING: Plotting',str(exfor_dict[index][2].shape[1])+'-column EXFOR data retrieved for subentry', exfor_dict[index][3]+', please make sure data look reasonable - column formatting is inconsistent for >4 columns.')
                        plt.errorbar(exfor_dict[index][2][:,4],exfor_dict[index][2][:,5], xerr=exfor_dict[index][2][:,0], yerr=exfor_dict[index][2][:,6], ls='none', capsize=1, label=exfor_dict[index][0]+' ('+exfor_dict[index][1]+')', marker=markers[k], markersize=4, linewidth=1, color=grey_colors[k])
                k=k+1

        except UnboundLocalError:
            print('Warning: Not able to plot exfor data!')

    elif exfor_manuel == True:
        path = f'../exfor/{foil_list[0][0:2]}/{reaction_product}_{xs_type}'
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
            label = f"{author_name} ({year})"

            plt.errorbar(energy, cross_section, xerr=energy_unc, yerr=cross_section_unc, ls='none', capsize=1, marker=markers[i], markersize=4, linewidth=1, color=grey_colors[i], label=label)



    if reaction_product == '86Y':
        plt.errorbar(E_86_haleema, xs_86_haleema, xerr=dE_86_haleema, yerr=d_xs_86_haleema, linewidth = 1, capsize = 1, marker='s', linestyle='None', color='black', label= 'Zaneb (2018)')
    if reaction_product == '87Y':
        plt.errorbar(E_87_haleema, xs_87_haleema, xerr=dE_87_haleema, yerr=d_xs_87_haleema, linewidth = 1, capsize = 1, marker='s', linestyle='None', color='black', label= 'Zaneb (2018)')
    if reaction_product == '88Y':
        plt.errorbar(E_88_haleema, xs_88_haleema, xerr=dE_88_haleema, yerr=d_xs_88_haleema, linewidth = 1, capsize = 1, marker='s', linestyle='None', color='black', label= 'Zaneb (2018)')


    plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='mediumseagreen', label='TALYS-2.0')
    plt.plot(energy_array, alice_spline(energy_array), linewidth=1, color='mediumblue', linestyle=(0,(5,5)), label='ALICE-2020')
    plt.plot(energy_array, coh_spline(energy_array), linewidth=1, color='olivedrab', linestyle=(0,(12,10)), label='CoH-3.5.3')
    plt.plot(energy_array, empire_spline(energy_array), linewidth=1.5, color='deepskyblue', linestyle=(0,(5,4,1,4,1,4)), label='EMPIRE 3.2.3')
    plt.plot(energy_array, tendl_spline(energy_array), linewidth=2, color='cornflowerblue', linestyle=':', label='TENDL-2023')
    
    plt.errorbar(df_88Zr_cum['Energy'], df_88Zr_cum['xs'], xerr=[df_88Zr_cum['Energy_min_unc'], df_88Zr_cum['Energy_plus_unc']], yerr=df_88Zr_cum['xs_unc'], marker='D', markersize=4, linewidth=2, linestyle='', color='hotpink', label='This work')

    plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
    plt.ylabel('Cross Section (mb)', fontsize=14)
    plt.title(f'{title}{xs_type}', fontsize=14)
    plt.legend(fontsize=12)
    plt.ylim((0))
    plt.xlim(0,50)
    if save_fig==True:
        plt.savefig(f'../Figures/PhD/xs_plots/{target_material}/{target_material}_dx_{reaction_product}_{state}_{xs_type}.pdf', dpi=600)
    plt.show()
    



#_________________________________Running the code______________________________________________________________________________________________________________

plot_xs('88ZR', 'ground_state', 40, 88, ['Zr08', 'Zr09', 'Zr10'], r'$^{nat}$Zr(d,x)$^{88}$Zr - ', 'Cumulative', save_fig=True, write_csv=False, exfor_manuel=False)
