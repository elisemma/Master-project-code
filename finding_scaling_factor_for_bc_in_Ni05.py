import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.optimize import curve_fit
from foil_class import Foil
import os 
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator


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



def calculate_chi2(observed, expected, unc_observed):
    observed = np.array(observed)
    expected = np.array(expected)
    diff = observed-expected 
    chi2 = np.sum(np.multiply(diff, diff)/np.multiply(unc_observed,unc_observed))
    return chi2



def run_chi2(data_by_reaction):
    observed_values = []
    expected_values = []
    unc_observed_values = []
    for reaction, data in data_by_reaction.items():
        #Getting monitor xs
        E_mon_list, xs_mon_list, xs_mon_unc_list = get_IAEA_monitro_xs(reaction)
        E_mon = np.array(E_mon_list)
        xs_mon = np.array(xs_mon_list)
        interp_xs_mon = PchipInterpolator(E_mon, xs_mon)

        #Getting observed data for the rection
        calc_xs = data['calc_xs']
        calc_xs_unc = data['calc_xs_unc']
        energies = data['energy']

        for i in range (len(calc_xs)):
            true_xs = interp_xs_mon(energies[i])
            observed_values.append(calc_xs[i])
            expected_values.append(true_xs)
            unc_observed_values.append(calc_xs_unc[i])
    dgf = 1 #This I have to think about

    chi2 = calculate_chi2(observed_values, expected_values, unc_observed_values)
    red_chi2 = chi2/(len(observed_values)-dgf) 
    return chi2, red_chi2















def caclulate_xs_in_stack(stack_df, compound, dp, scaling_factor):
    # stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.2f}.csv')
   

    monitor_stack_df = stack_df[stack_df['name']=='Ni05']


    calc_xs_list_of_list = [] #This list will cointain lists of cross sections calculated from beam currents for all the mon reactions in the foils: [[Ni01:56CO, Ni01:58CO, Ni01:61CU], [Ni02:56CO, Ni02:58CO, Ni02:61CU], ...]
    calc_xs_unc_list_of_list = [] #Same as calc_xs_list_of_list but with the uncertainties

    beam_energy_in_foil_list_list = []

    reaction_list_list = []


    for index, row in monitor_stack_df.iterrows():
        foil_name = row['name']
        target_material = row['compound']
        # beam_energy_in_foil = row['mu_E']

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

        foil = Foil(foil_name, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, dp, scaling_factor)

        foil.assign_molar_mass()
        foil.assign_areal_dens_unc_percent() 
        foil.calculate_decay_constant()
        foil.find_monitor_cross_section()
        foil.calculate_beam_currents_w_unc()
        foil.calculate_weighted_average_beam_current()
        foil.convert_beam_current_back_to_xs_w_unc()

        beam_energy_in_foil = foil.beam_energy_in_foil
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












def plot_chi2(scaling_factor_list):
    chi2_list = []
    red_chi2_list = []
    dp = 0.977


    for scaling_factor in scaling_factor_list:
        # calculate beam current here
        beam_current_list=[]
        beam_current_unc_list=[]
        energy_list=[]

        # for foil in compartment:
        stack_df = pd.read_csv(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{dp:.3f}.csv')

        monitor_compounds = ['Ni']
        compartment_list = ['05']
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

           

            data_by_reaction = {}

            for compound in monitor_compounds:
                data_dict = caclulate_xs_in_stack(monitor_stack_df, compound, dp, scaling_factor)
                data_by_reaction.update(data_dict)
               
        chi2, red_chi2 = run_chi2(data_by_reaction)
        chi2_list.append(chi2)
        red_chi2_list.append(red_chi2)

    min_index = red_chi2_list.index(min(red_chi2_list))
    min_scaling_factor = scaling_factor_list[min_index]


    plt.plot(scaling_factor_list, red_chi2_list)
    plt.xlabel('scaling factor')
    plt.ylabel('reduced chi2')
    plt.title(f'Minimized when scaling factor = {min_scaling_factor:.3f}')
    plt.show()
    print(monitor_stack_df)




scaling_factor_list = np.linspace(1,90,10000)
plot_chi2(scaling_factor_list)








# data_by_reaction = {'58CO': {'beam_current': [38.12635144593368], 'beam_current_unc': [8.854036377386427], 'energy': [2.2266116583208797]},
                    # '61CU': {'beam_current': [0.5806523386626631], 'beam_current_unc': [0.13051120241881975], 'energy': [2.2266116583208797]}}