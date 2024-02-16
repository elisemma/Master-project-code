
import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.optimize import curve_fit
from foil_class import Foil
import os 
from scipy.interpolate import interp1d



def get_IAEA_monitro_xs(reaction_product):
    filename = './Monitor_cross_section_data/IAEA_monitor_xs_' + reaction_product + '.txt'
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

    return E_mon_list, xs_list, xs_unc_list


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

    elif method == 'global_xs':
        print('xs')
        data_by_reaction = x_data

        y_data = []
        true_array = []
        unc_data = []

        for reaction, data in data_by_reaction.items():
            #Getting monitor xs
            E_mon_list, xs_mon_list, xs_mon_unc_list = get_IAEA_monitro_xs(reaction)
            E_mon = np.array(E_mon_list)
            xs_mon = np.array(xs_mon_list)
            interp_xs_mon = interp1d(E_mon, xs_mon,kind='linear')

            #Getting observed data for the rection
            calc_xs = data['calc_xs']
            calc_xs_unc = data['calc_xs_unc']
            energies = data['energy']

            for i in range (len(calc_xs)):
                true_xs = interp_xs_mon(energies[i])
                y_data.append(calc_xs[i])
                true_array.append(true_xs)
                unc_data.append(calc_xs_unc[i])
        dgf = 1 #This I have to think about

    else:
        print(f'ERROR: The {method} method for variance minimization is not implemented.')

    chi2 = calculate_chi2(y_data, true_array, unc_data)
    red_chi2 = chi2/(len(y_data)-dgf) 
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


def caclulate_xs_in_stack(stack_df, compound, dp):
    # stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.2f}.csv')
   

    monitor_stack_df = stack_df[stack_df['name'].str.contains(f'{compound}')]

    calc_xs_list_of_list = [] #This list will cointain lists of cross sections calculated from beam currents for all the mon reactions in the foils: [[Ni01:56CO, Ni01:58CO, Ni01:61CU], [Ni02:56CO, Ni02:58CO, Ni02:61CU], ...]
    calc_xs_unc_list_of_list = [] #Same as calc_xs_list_of_list but with the uncertainties

    beam_energy_in_foil_list_list = []

    reaction_list_list = []


    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']
        beam_energy_in_foil = row['mu_E']

        A0_by_curie_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_curie.csv')

        if os.path.exists(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv'):
            A0_by_hand_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Calculated_A0/{foil_name}_A0_by_hand.csv')
            A0_concat_df = pd.concat((A0_by_hand_df, A0_by_curie_df), axis=0)

        else:
            A0_concat_df = A0_by_curie_df


        # print('A0_conccat: ', A0_concat_df)

        reaction_list = A0_concat_df['Isotope'].tolist()
        A0_list = A0_concat_df['A0'].tolist()
        A0_unc_list = A0_concat_df['A0_unc'].tolist()
        A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]
        areal_dens = row['areal_density']
        areal_dens_unc_percent = 0.2 #XXXXXXXXXXXX this is not true, need to find it

        foil = Foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent, dp)

        foil.assign_molar_mass()
        foil.calculate_decay_constant()
        foil.find_monitor_cross_section()
        foil.calculate_beam_currents_w_unc()
        foil.calculate_weighted_average_beam_current()
        foil.convert_beam_current_back_to_xs_w_unc()

        foil_calc_xs_list = foil.calc_xs_list
        foil_calc_xs_unc_list = foil.calc_xs_unc_list


        beam_energy_in_foil_array = np.zeros(len(reaction_list))
        beam_energy_in_foil_array.fill(beam_energy_in_foil)
        beam_energy_in_foil_list_list.append(beam_energy_in_foil_array)
        reaction_list_list.append(reaction_list)
        calc_xs_list_of_list.append(foil_calc_xs_list)
        calc_xs_unc_list_of_list.append(foil_calc_xs_unc_list)



        # Initialize nested dictionary to store beam currents, uncertainties and energies for each reaction
        data_by_reaction = {}

        # Iterate through reaction_list_list
        for i, reaction_list in enumerate(reaction_list_list):
            for j, reaction in enumerate(reaction_list):
                if reaction not in data_by_reaction:
                    data_by_reaction[reaction] = {'calc_xs': [], 'calc_xs_unc': [], 'energy': []}  # Initialize inner dictionary with lists
                # Add calc_xs, calc_xs_unc and energy corresponding to the reaction
                data_by_reaction[reaction]['calc_xs'].append(calc_xs_list_of_list[i][j])
                data_by_reaction[reaction]['calc_xs_unc'].append(calc_xs_unc_list_of_list[i][j])
                data_by_reaction[reaction]['energy'].append(beam_energy_in_foil_list_list[i][j])

    return data_by_reaction



def plot_chi2(dp_list, compartment_list, method):
    chi2_list = []
    red_chi2_list = []
    method_list = ['p0', 'p1']


    for dp in dp_list:
        # calculate beam current here
        print(f'dp={dp}_________________________')
        beam_current_list=[]
        beam_current_unc_list=[]
        energy_list=[]

        # for foil in compartment:
        stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.3f}.csv')

        monitor_compounds = ['Ni', 'Ti']
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

            A0_list = A0_concat_df['A0'].tolist()
            A0_unc_list = A0_concat_df['A0_unc'].tolist()
            A0_unc_list = [1e10 if np.isinf(value) else value for value in A0_unc_list]


            if method in method_list:

                foil_beam_cur_list, foil_beam_cur_unc_list = beam_currents_in_foil(foil_name, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent, dp)

                for beam_current, beam_current_unc in zip(foil_beam_cur_list, foil_beam_cur_unc_list):
                    beam_current_list.append(beam_current)
                    beam_current_unc_list.append(beam_current_unc)
                    energy_list.append(beam_energy_in_foil)
                x_data = energy_list
                y_data = beam_current_list
                unc_data = beam_current_unc_list


            elif method == 'global_xs':

                data_by_reaction = {}

                for compound in monitor_compounds:
                    data_dict = caclulate_xs_in_stack(monitor_stack_df, compound, dp)
                    data_by_reaction.update(data_dict)
                x_data = data_by_reaction
                y_data = []
                unc_data = []
               
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




#_____________________Running the code___________________________

dp_array1 = np.arange(0.8, 1.21, 0.01)
dp_array2 = np.arange(0.97, 0.99, 0.001)
dp_array3 = np.union1d(dp_array1, dp_array2)
dp_array4 = np.arange(0.99, 1.001, 0.001)
dp_array5 = np.union1d(dp_array3, dp_array4)
dp_array6 = np.arange(0.94, 0.96, 0.001)
dp_array = np.union1d(dp_array5, dp_array6)

# plot_chi2(dp_array, ['01', '02', '03', '04'], 'global_xs')
plot_chi2(dp_array, ['04'], 'p0')


print(dp_array)




















