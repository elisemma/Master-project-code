
import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.optimize import curve_fit
from foil_class import Foil
import os 


# Plan: Get the beam currents with uncertainty for one compartment
# Then make a p0 fit whith these data points
# Give the beam currents as observed and the p0 fit as the true value to chi2 
# plot the chi2 as a function of dp and find the minimum 

def p0(x, a):
    return a


def fit_p0(x_data, y_data, unc_data):
    y_data = np.array(y_data)
    y_data[np.isnan(y_data)] = 0
    popt, cov  = curve_fit(p0, x_data, y_data, p0=100, sigma=unc_data)
    return popt[0]


def calculate_chi2(observed, expected, unc_observed):
    diff = observed-expected 
    chi2 = np.sum(np.multiply(diff, diff)/np.multiply(unc_observed,unc_observed))
    return chi2


def run_chi2(x_data, y_data, unc_data):
    true = fit_p0(x_data, y_data, unc_data)
    true_array = np.zeros(len(y_data))
    true_array.fill(true)
    chi2 = calculate_chi2(y_data, true_array, unc_data)
    red_chi2 = chi2/(len(y_data)-1) # NB: Only -1 when using p0
    print('true: ', true,', y_data[0]: ', y_data[0])

    # plt.errorbar(x_data, y_data, yerr=unc_data, color='hotpink', label='data', marker='*', linestyle='None')
    # plt.plot(x_data, true_array, color='skyblue', label='curve fit')
    # plt.text(0.05, 0.95, f'chi2 = {chi2}', transform=plt.gca().transAxes, fontsize=12,verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    # plt.legend()
    # plt.show()
    return chi2, red_chi2


def beam_currents_in_foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent, dp):
    foil = Foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent, dp)

    foil.assign_molar_mass()
    foil.calculate_decay_constant()
    foil.find_monitor_cross_section()
    foil.calculate_beam_currents_w_unc()

    foil_beam_cur_list = foil.beam_current_list
    foil_beam_cur_unc_list = foil.beam_current_unc_list

    return foil_beam_cur_list, foil_beam_cur_unc_list


def plot_chi2(dp_list, compartment):
    chi2_list = []
    red_chi2_list = []

    for dp in dp_list:
        # calculate beam current here
        print(f'dp={dp}_________________________')
        beam_current_list=[]
        beam_current_unc_list=[]
        energy_list=[]

        # for foil in compartment:
        stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.2f}.csv')

        monitor_compounds = ['Ni', 'Ti']
        monitor_stack_df = pd.concat([stack_df[(stack_df['name'].str.contains(compound)) & (stack_df['name'].str.contains(compartment))] for compound in monitor_compounds])


        print(monitor_stack_df)

        for index, row in monitor_stack_df.iterrows():
            foil_name = row['name']
            target_material = row['compound']
            beam_energy_in_foil = row['mu_E']
            areal_dens = row['areal_density']
            areal_dens_unc_percent = 2 #XXXXXXXXXXXX this is not true, need to find it

            A0_by_curie_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_curie.csv')

            if os.path.exists(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv'):
                A0_by_hand_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv')
                A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)


            else:
                A0_concat_df = A0_by_curie_df

            reaction_list = A0_concat_df['Isotope'].tolist()
            # if target_material =='Ni':
            #     reaction_list = ['58CO']
            # if target_material == 'Ti':
            #     reaction_list = ['48V']

            # A0_concat_df = A0_concat_df[A0_concat_df['Isotope']==reaction_list[0]]

            A0_list = A0_concat_df['A0'].tolist()
            A0_unc_list = A0_concat_df['A0_unc'].tolist()
            A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]


            foil_beam_cur_list, foil_beam_cur_unc_list = beam_currents_in_foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent, dp)

            for beam_current, beam_current_unc in zip(foil_beam_cur_list, foil_beam_cur_unc_list):
                beam_current_list.append(beam_current)
                beam_current_unc_list.append(beam_current_unc)
                energy_list.append(beam_energy_in_foil)

        chi2, red_chi2 = run_chi2(energy_list, beam_current_list,beam_current_unc_list)
        chi2_list.append(chi2)
        red_chi2_list.append(red_chi2)

    min_index = red_chi2_list.index(min(red_chi2_list))
    min_dp = dp_list[min_index]

    plt.plot(dp_list, red_chi2_list)
    plt.xlabel('dp')
    plt.ylabel('reduced chi2')
    plt.title(f'Compartment {compartment}, minimized when dp = {min_dp:.2f}')
    plt.show()


dp_array = np.arange(0.8, 1.21, 0.01)
# dp_array=[0.90, 1.00, 1.10]

plot_chi2(dp_array, '05')




















