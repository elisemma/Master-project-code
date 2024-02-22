import numpy as np  
import matplotlib.pyplot as plt 
import pandas as pd 
from foil_class_for_xs_calc import Foil 
from scipy.interpolate import splrep, splev
import os 
import csv



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


def get_talys_data(filename):
    E_talys, xs_talys = np.loadtxt('./talys_calculations/'+filename, unpack = True, skiprows = 24)
    E_list = [0]+list(E_talys)
    xs_list = [0]+list(xs_talys)
    return E_list, xs_list



def plot_xs(reaction_product, Z, A, foil_list): #NB this function assumes that the first foil in the foil list is Zr01 and that you don't skip any foils. 

    Zr_energy_data = {'Zr01': {'energy': 26.375581480000005, 'min_unc': 0.605581480000005, 'plus_unc': 0.5944185199999943},
                    'Zr02': {'energy': 19.1745053, 'min_unc': 0.7845052999999993, 'plus_unc': 0.7754946999999994}, 
                    'Zr03': {'energy': 12.309193775999999, 'min_unc': 1.059193775999999, 'plus_unc': 1.1008062240000012}, 
                    'Zr04': {'energy': 6.082350534817555, 'min_unc': 1.4923505348175548, 'plus_unc': 1.8676494651824447}, 
                    'Zr05': {'energy': 1.5675066930508672, 'min_unc': 1.2375066930508674, 'plus_unc': 0.7424933069491328}}

    Zr_energy_list = [Zr_energy_data[foil]['energy'] for foil in foil_list]
    Zr_energy_min_unc_list = [Zr_energy_data[foil]['min_unc'] for foil in foil_list]
    Zr_energy_plus_unc_list = [Zr_energy_data[foil]['plus_unc'] for foil in foil_list]

    xs_list = []
    xs_unc_list = []

    for foil in foil_list:
        xs, xs_unc = calc_xs(foil, reaction_product)
        xs_list.append(xs)
        xs_unc_list.append(xs_unc)

    talys_file = f'rp0{Z}0{A}.tot'
    E_talys, xs_talys = get_talys_data(talys_file)
    talys_spline = splrep(E_talys, xs_talys)
    energy_array = np.linspace(0,30,300)


    csv_file_path = f'./Calculated_xs/Zr_foils/{reaction_product}_xs.csv'
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Foil', 'Energy', 'Energy_min_unc', 'Energy_plus_unc', 'xs', 'xs_unc'])
        for i in range(len(foil_list)):
            csv_writer.writerow([f'{foil_list[i]}', f'{Zr_energy_list[i]}', f'{Zr_energy_min_unc_list[i]}', f'{Zr_energy_plus_unc_list[i]}', f'{xs_list[i]}', f'{xs_unc_list[i]}'])

    plt.plot(energy_array, splev(energy_array, talys_spline), color='gold', label='TALYS')

    plt.errorbar(Zr_energy_list[0:len(foil_list)], xs_list, xerr=[Zr_energy_min_unc_list[0:len(foil_list)], Zr_energy_plus_unc_list[0:len(foil_list)]], yerr=xs_unc_list, marker='d', markersize=5, linestyle='', color='hotpink', label='This work')
    plt.xlabel('Beam energy (MeV)')
    plt.ylabel('Cross section (mb)')
    plt.title(reaction_product)
    plt.legend()
    plt.xlim(0,30)
    plt.show()





plot_xs('96NB', 41, 96, ['Zr01', 'Zr02', 'Zr03', 'Zr04'])
plot_xs('90NB', 41, 90, ['Zr01', 'Zr02', 'Zr03', 'Zr04'])





