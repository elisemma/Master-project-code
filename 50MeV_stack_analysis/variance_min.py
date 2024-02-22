
import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.optimize import curve_fit
from foil_class import Foil
import os 
from scipy.interpolate import interp1d




def p0(x, a):
    return a


def p1(x, a, b):
    x = np.array(x)
    b_array = np.zeros(len(x))
    b_array.fill(b)
    return a*x + b_array


def fit_p0(x_data, y_data, unc_data):
    y_data = np.array(y_data)
    y_data[np.isnan(y_data)] = 0
    popt, cov  = curve_fit(p0, x_data, y_data, p0=100, sigma=unc_data)
    return popt[0]


def fit_p1(x_data, y_data, unc_data):
    y_data = np.array(y_data)
    y_data[np.isnan(y_data)] = 0
    popt, cov  = curve_fit(p1, x_data, y_data, p0=[0, 100], sigma=unc_data)
    return popt[0], popt[1]


def calculate_chi2(observed, expected, unc_observed):
    observed = np.array(observed)
    expected = np.array(expected)
    diff = observed-expected 
    chi2 = np.sum(np.multiply(diff, diff)/np.multiply(unc_observed,unc_observed))
    return chi2


def run_chi2(x_data, y_data, unc_data, method):

    if method == 'p0':
        print('p0')
        true = fit_p0(x_data, y_data, unc_data)
        true_array = np.zeros(len(y_data))
        true_array.fill(true)
        dgf = 1

    elif method == 'p1':
        print('p1')
        a,b = fit_p1(x_data, y_data, unc_data)
        x = np.array(x_data)
        b_array = np.zeros(len(x))
        b_array.fill(b)
        true_array = a*x+b_array
        dgf = 2

    else:
        print(f'ERROR: The {method} method for variance minimization is not implemented.')

    chi2 = calculate_chi2(y_data, true_array, unc_data)
    red_chi2 = chi2/(len(y_data)-dgf) 
    return chi2, red_chi2


def beam_currents_in_foil(foil_name, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, dp):
    foil = Foil(foil_name, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, dp)

    foil.assign_molar_mass()
    foil.assign_areal_dens_unc_percent()
    foil.calculate_decay_constant()
    foil.find_monitor_cross_section()
    foil.calculate_beam_currents_w_unc()

    foil_beam_cur_list = foil.beam_current_list
    foil_beam_cur_unc_list = foil.beam_current_unc_list
    beam_energy_in_foil = foil.beam_energy_in_foil

    return beam_energy_in_foil, foil_beam_cur_list, foil_beam_cur_unc_list


def plot_chi2(dp_list, compartment_list, method):
    chi2_list = []
    red_chi2_list = []

    for dp in dp_list:
        # calculate beam current here
        beam_current_list=[]
        beam_current_unc_list=[]
        energy_list=[]

        # for foil in compartment:
        stack_df = pd.read_csv(f'./Stack_calculations/stack_50MeV_dp_{dp:.3f}.csv')
        monitor_compounds = ['Fe', 'Ti']
        # if compartment == '05':
        #     monitor_compounds = ['Ni']
        # monitor_stack_df = pd.concat([stack_df[(stack_df['name'].str.contains(compound)) & (stack_df['name'].str.contains(compartment))] for compound in monitor_compounds])
        monitor_stack_df = pd.concat([stack_df[(stack_df['name'].str.contains(compound)) & 
                                      (stack_df['name'].str.contains(compartment))] 
                              for compound in monitor_compounds 
                              for compartment in compartment_list])

         
        # print(monitor_stack_df)

        for index, row in monitor_stack_df.iterrows():
            foil_name = row['name']
            target_material = row['compound']
            # beam_energy_in_foil = row['mu_E']
            areal_dens = row['areal_density']

            A0_by_curie_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_curie.csv')

            A0_concat_df = A0_by_curie_df

            reaction_list = A0_concat_df['Isotope'].tolist()

            A0_list = A0_concat_df['A0'].tolist()
            A0_unc_list = A0_concat_df['A0_unc'].tolist()
            A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]


            

            beam_energy_in_foil, foil_beam_cur_list, foil_beam_cur_unc_list = beam_currents_in_foil(foil_name, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, dp)
            if np.isnan(beam_energy_in_foil):
                beam_energy_in_foil = 0

            for beam_current, beam_current_unc in zip(foil_beam_cur_list, foil_beam_cur_unc_list):
                beam_current_list.append(beam_current)
                beam_current_unc_list.append(beam_current_unc)
                energy_list.append(beam_energy_in_foil)
            x_data = energy_list
            y_data = beam_current_list
            unc_data = beam_current_unc_list

        chi2, red_chi2 = run_chi2(x_data, y_data, unc_data, method)
        chi2_list.append(chi2)
        red_chi2_list.append(red_chi2)

    min_index = red_chi2_list.index(min(red_chi2_list))
    min_dp = dp_list[min_index]


    plt.plot(dp_list, red_chi2_list)
    plt.xlabel('dp')
    plt.ylabel('reduced chi2')
    plt.title(f'Compartment {compartment_list}, method: {method}, minimized when dp = {min_dp:.3f}')
    plt.show()
    print(monitor_stack_df)




#_____________________Running the code___________________________

dp_array1 = np.arange(0.8, 1.21, 0.01)
dp_array2 = np.arange(0.95, 0.98, 0.001)
dp_array = np.union1d(dp_array1, dp_array2)

plot_chi2(dp_array, ['05', '11'], 'p1')


# print(dp_array)
