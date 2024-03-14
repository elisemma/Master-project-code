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





def plot_xs(reaction_product, state, Z, A, foil_list, title, write_csv=False, save_fig=False, use_R_to_calc_xs=True):

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
            elif state == 'isomeric_state':
                xs_state = 'm1'

            reaction_product_w_state = reaction_product+xs_state

            xs, xs_unc = calc_xs_from_R(foil, reaction_product, reaction_product_w_state)
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
        coh_file = f'0{Z}-0{A}{reaction_product[2:]}{coh_state}_coh.txt'
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



    #Get data from exfor
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

    except UnboundLocalError:
        print('Warning: Not able to plot exfor data!')


    

    # plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='gold', label='TALYS')
    # plt.plot(energy_array, splev(energy_array, alice_spline), linewidth=1, color='plum', linestyle='--', label='ALICE')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.2, linewidth=1, color='skyblue', label='Test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*0.9, linewidth=1, color='mediumaquamarine', linestyle='-.', label='test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.1, linewidth=2, color='mediumpurple', linestyle=':', label='test.')

    plt.plot(energy_array, splev(energy_array, talys_spline), linewidth=1, color='mediumaquamarine', label='TALYS-2.0')
    plt.plot(energy_array, alice_spline(energy_array), linewidth=1, color='deepskyblue', linestyle='--', label='ALICE-2020')
    plt.plot(energy_array, coh_spline(energy_array), linewidth=1, color='skyblue', label='CoH-3.5.3')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*0.9, linewidth=1, color='darkturquoise', linestyle='-.', label='test')
    # plt.plot(energy_array, splev(energy_array, talys_spline)*1.1, linewidth=2, color='cornflowerblue', linestyle=':', label='test.')


    # plt.errorbar(energy_list[0:len(foil_list)], xs_list, xerr=[energy_min_unc_list[0:len(foil_list)], energy_plus_unc_list[0:len(foil_list)]], yerr=xs_unc_list, marker='D', markersize=5, linewidth=3, linestyle='', color='hotpink', label='This work')
    plt.errorbar(energy_list, xs_list, xerr=[energy_min_unc_list, energy_plus_unc_list], yerr=xs_unc_list, marker='D', markersize=5, linewidth=3, linestyle='', color='hotpink', label='This work')

    plt.xlabel('Deuteron Energy (MeV)', fontsize=14)
    plt.ylabel('Cross Section (mb)', fontsize=14)
    plt.title(title, fontsize=14)
    plt.legend(fontsize=14)
    plt.ylim((0))
    plt.xlim(0,50)
    if save_fig==True:
        plt.savefig(f'./Figures/xs_plots/{target_material}_dx_{reaction_product}.pdf', dpi=600)
    plt.show()
    



#_________________________________Running the code_______________________________________


#_______________________________Zr-foils________________________________________________
plot_xs('86Y', 'ground_state', 39, 86, ['Zr06', 'Zr07', 'Zr08', 'Zr09'], r'$^{nat}$Zr(d,x)$^{86}$Y - Cumulative', save_fig=True)

# plot_xs('87Y', '87Yg', 39, 87, ['Zr01', 'Zr02', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('87Ym', '87Ym1', 39, 87, ['Zr01', 'Zr02', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])

# plot_xs('88Y', '88Yg', 39, 88, ['Zr01', 'Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('88ZR', '88ZRg', 40, 88, ['Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('88NB', '88NBg', 41, 88, ['Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('88Y', '88Yg', 39, 88, ['Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])

# plot_xs('89NB', '89NBg', 41, 89, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('89ZR', '89ZRg', 40, 89, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('89Ym', '89Ym1', 39, 89, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])

# plot_xs('90NB', '90NBg', 41, 90, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])

# plot_xs('90Ym', '90Ym1', 39, 90, ['Zr01', 'Zr03', 'Zr06', 'Zr07', 'Zr08', 'Zr10'])

# plot_xs('91Y', '91Yg', 39, 91, ['Zr03', 'Zr04', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])

# plot_xs('92Y', '92Yg', 39, 92, ['Zr01', 'Zr03', 'Zr06', 'Zr07', 'Zr08'])

# plot_xs('95ZR', '95ZRg', 40, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('95NBm', '95NBm1', 41, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])
# plot_xs('95NB', '95NBg', 41, 95, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])

# plot_xs('96NB', '96NBg', 41, 96, ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05', 'Zr06', 'Zr07', 'Zr08', 'Zr09', 'Zr10'])



#_______________________________Ni-foils________________________________________________

# plot_xs('52MN', '52MNg', 25, 52, ['Ni01', 'Ni02'])

# plot_xs('54MN', '54MNg', 25, 54, ['Ni01'])

# plot_xs('55CO', '55COg', 27, 55, ['Ni01', 'Ni02', 'Ni03', 'Ni05'])

# plot_xs('56CO', '56COg', 27, 56, ['Ni01', 'Ni02', 'Ni03', 'Ni04'])

# plot_xs('57NI', '57NIg', 28, 57, ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05'])
# plot_xs('57CO', '57COg', 27, 57, ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05'])

# plot_xs('58CO', '58COg', 27, 58, ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05'])

# plot_xs('60CO', '60COg', 27, 60, ['Ni01', 'Ni04'])

# plot_xs('60CU', '60CUg', 29, 60, ['Ni02', 'Ni03'])

# plot_xs('61CU', '61CUg', 29, 61, ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05'])

# plot_xs('64CU', '64CUg', 29, 64, ['Ni02', 'Ni03', 'Ni04'])

# plot_xs('65NI', '65NIg', 28, 65, ['Ni01', 'Ni02', 'Ni03', 'Ni04'])





#_______________________________Ti-foils________________________________________________

# plot_xs('46SC', '46SCg', 21, 46, ['Ti01', 'Ti02', 'Ti03', 'Ti04', 'Ti05', 'Ti06', 'Ti08', 'Ti09', 'Ti10', 'Ti11'])

# plot_xs('47SC', '47SCg', 21, 47, ['Ti01', 'Ti02', 'Ti03', 'Ti06', 'Ti08', 'Ti09', 'Ti10', 'Ti11'])

# plot_xs('48SC', '48SCg', 21, 48, ['Ti06', 'Ti08', 'Ti09', 'Ti10'])

# plot_xs('48V', '48Vg', 23, 48, ['Ti01', 'Ti02', 'Ti03', 'Ti04', 'Ti05', 'Ti06', 'Ti08', 'Ti09', 'Ti10', 'Ti11'])



















