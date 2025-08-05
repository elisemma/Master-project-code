
import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import numpy as np 



# Fe peak data _____________________________________________________________________________________________________________________
path_Fe = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/50MeV_stack_analysis/MyGeneratedFiles/Fe_foils/'
# 50 MeV:
file_BUFE01 = 'BU170217_Fe01_18cm_50MeV/BU170217_Fe01_18cm_50MeV_peak_data.csv'
file_BVFE02 = 'BV170217_Fe02_18cm_50MeV/BV170217_Fe02_18cm_50MeV_peak_data.csv'
file_BWFE03 = 'BW170217_Fe03_18cm_50MeV/BW170217_Fe03_18cm_50MeV_peak_data.csv'
file_BYFE04 = 'BY190217_Fe04_18cm_50MeV/BY190217_Fe04_18cm_50MeV_peak_data.csv'
file_BZFE05 = 'BZ190217_Fe05_18cm_50MeV/BZ190217_Fe05_18cm_50MeV_peak_data.csv'


file_concat_Fe = 'combined_peak_data_Fe.csv'


df_BUFE01 = pd.read_csv(path_Fe+file_BUFE01)
df_BVFE02 = pd.read_csv(path_Fe+file_BVFE02)
df_BWFE03 = pd.read_csv(path_Fe+file_BWFE03)
df_BYFE04 = pd.read_csv(path_Fe+file_BYFE04)
df_BZFE05 = pd.read_csv(path_Fe+file_BZFE05)







df_concat_Fe = pd.concat((df_BUFE01, df_BVFE02, df_BWFE03, df_BYFE04, df_BZFE05), axis = 0)
df_concat_Fe.to_csv(path_Fe+file_concat_Fe)





pd.set_option('display.max_rows', None)
df_57CO=df_concat_Fe[df_concat_Fe['isotope']=='57CO']



# Get unique energy values and sort them
sorted_unique_energies = df_57CO['energy'].drop_duplicates().sort_values()
# Optionally convert to a list
sorted_unique_energies_list = sorted_unique_energies.tolist()
# Print them
print(sorted_unique_energies_list)








def fit_prod_rate(isotope_list, isotope_chain_parent, foil, path, file, stack, plot=False):
    t_irr_h = 0.33

    R = {}
    for isotope in isotope_list:
        if isotope == '86ZR':
            R0 = 0
        else:
            R0 = 1E4
        R[isotope] =[[R0,t_irr_h]] 

    if stack == '30MeV':
        EoB = '02/13/2017 14:47:00'

    if stack == '50MeV':
        EoB = '02/12/2017 19:21:00'

    dc = ci.DecayChain(isotope_chain_parent, R=R, units='h')
    dc.get_counts(foil, EoB, path+file)
    isotopes, R, cov_R = dc.fit_R()
    print(f'A0 for {isotope_chain_parent}: ', dc.activity(isotope_chain_parent, 0))

    if plot==True:
        dc.plot()

    return isotopes, R, cov_R




def calc_prod_rates_in_foil(isotope_list_list, isotope_chain_parent_list, foil, stack, write_to_file=False, show_plot=False):
    print(' ')
    print(foil)
    csv_file_path = f'./Calculated_R/{foil}_R_by_curie.csv'
    if write_to_file==True:
            with open(csv_file_path, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerow(['Isotope', 'R', 'R_unc', 'Rel_R_unc_in_percent'])

    for isotope_list, parent in zip(isotope_list_list, isotope_chain_parent_list):
        isotopes, R, cov_R = fit_prod_rate(isotope_list, parent, foil, path_Fe, file_concat_Fe, stack, plot=show_plot)

        print(isotopes)
        print(R)
        # print(cov_R)
        if write_to_file==True:
            with open(csv_file_path, 'a', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                for i, isotope in enumerate(isotopes):
                    # print(i, isotope)
                    csv_writer.writerow([f'{isotope}', f'{R[i]}', f'{np.sqrt(cov_R[i,i])}', f'{np.sqrt(cov_R[i,i])/R[i]*100:.2e}'])




# isotope_list_list_Fe01 = [['48V']]
# isotope_chain_parent_list_Fe01 = ['48V']

# isotope_list_list_Fe02 = [['48V']]
# isotope_chain_parent_list_Fe02 = ['48V']

# isotope_list_list_Fe03 = [['48V']]
# isotope_chain_parent_list_Fe03 = ['48V']

# isotope_list_list_Fe04 = [['48V']]
# isotope_chain_parent_list_Fe04 = ['48V']

# isotope_list_list_Fe05 = [['48V']]
# isotope_chain_parent_list_Fe05 = ['48V']
       

# isotope_list_list_Fe01 = [['52MN']]
# isotope_chain_parent_list_Fe01 = ['52MN']

# isotope_list_list_Fe02 = [['52MN']]
# isotope_chain_parent_list_Fe02 = ['52MN']

# isotope_list_list_Fe03 = [['52MN']]
# isotope_chain_parent_list_Fe03 = ['52MN']

# isotope_list_list_Fe04 = [['52MN']]
# isotope_chain_parent_list_Fe04 = ['52MN']

# isotope_list_list_Fe05 = [['52MN']]
# isotope_chain_parent_list_Fe05 = ['52MN']



# isotope_list_list_Fe01 = [['54MN']]
# isotope_chain_parent_list_Fe01 = ['54MN']

# isotope_list_list_Fe02 = [['54MN']]
# isotope_chain_parent_list_Fe02 = ['54MN']

# isotope_list_list_Fe03 = [['54MN']]
# isotope_chain_parent_list_Fe03 = ['54MN']

# isotope_list_list_Fe04 = [['54MN']]
# isotope_chain_parent_list_Fe04 = ['54MN']

# isotope_list_list_Fe05 = [['54MN']]
# isotope_chain_parent_list_Fe05 = ['54MN']



# isotope_list_list_Fe01 = [['55CO']]
# isotope_chain_parent_list_Fe01 = ['55CO']

# isotope_list_list_Fe02 = [['55CO']]
# isotope_chain_parent_list_Fe02 = ['55CO']

# isotope_list_list_Fe03 = [['55CO']]
# isotope_chain_parent_list_Fe03 = ['55CO']

# isotope_list_list_Fe04 = [['55CO']]
# isotope_chain_parent_list_Fe04 = ['55CO']

# isotope_list_list_Fe05 = [['55CO']]
# isotope_chain_parent_list_Fe05 = ['55CO']
       


# isotope_list_list_Fe01 = [['56CO']]
# isotope_chain_parent_list_Fe01 = ['56CO']

# isotope_list_list_Fe02 = [['56CO']]
# isotope_chain_parent_list_Fe02 = ['56CO']

# isotope_list_list_Fe03 = [['56CO']]
# isotope_chain_parent_list_Fe03 = ['56CO']

# isotope_list_list_Fe04 = [['56CO']]
# isotope_chain_parent_list_Fe04 = ['56CO']

# isotope_list_list_Fe05 = [['56CO']]
# isotope_chain_parent_list_Fe05 = ['56CO']



# isotope_list_list_Fe01 = [['57CO']]
# isotope_chain_parent_list_Fe01 = ['57CO']

# isotope_list_list_Fe02 = [['57CO']]
# isotope_chain_parent_list_Fe02 = ['57CO']

# isotope_list_list_Fe03 = [['57CO']]
# isotope_chain_parent_list_Fe03 = ['57CO']

# isotope_list_list_Fe04 = [['57CO']]
# isotope_chain_parent_list_Fe04 = ['57CO']

# isotope_list_list_Fe05 = [['57CO']]
# isotope_chain_parent_list_Fe05 = ['57CO']




# isotope_list_list_Fe03 = [['58CO']]
# isotope_chain_parent_list_Fe03 = ['58CO']

# isotope_list_list_Fe04 = [['58CO']]
# isotope_chain_parent_list_Fe04 = ['58CO']

# isotope_list_list_Fe05 = [['58CO']]
# isotope_chain_parent_list_Fe05 = ['58CO']










# isotope_list_list_Fe01 = [['48V'], ['52MN'], ['54MN'], ['55CO'], ['56CO'], ['57CO']]
# isotope_chain_parent_list_Fe01 = ['48V', '52MN', '54MN', '55CO', '56CO', '57CO']

# isotope_list_list_Fe02 = [['48V'], ['52MN'], ['54MN'], ['55CO'], ['56CO'], ['57CO']]
# isotope_chain_parent_list_Fe02 = ['48V', '52MN', '54MN', '55CO', '56CO', '57CO']

# isotope_list_list_Fe03 = [['48V'], ['52MN'], ['54MN'], ['55CO'], ['56CO'], ['57CO'], ['58CO']]
# isotope_chain_parent_list_Fe03 = ['48V', '52MN', '54MN', '55CO', '56CO', '57CO', '58CO']

# isotope_list_list_Fe04 = [['48V'], ['52MN'], ['54MN'], ['55CO'], ['56CO'], ['57CO'], ['58CO']]
# isotope_chain_parent_list_Fe04 = ['48V', '52MN', '54MN', '55CO', '56CO', '57CO', '58CO']

# isotope_list_list_Fe05 = [['48V'], ['52MN'], ['54MN'], ['55CO'], ['56CO'], ['57CO'], ['58CO']]
# isotope_chain_parent_list_Fe05 = ['48V', '52MN', '54MN', '55CO', '56CO', '57CO', '58CO']




# calc_prod_rates_in_foil(isotope_list_list_Fe01, isotope_chain_parent_list_Fe01, 'Fe01', '50MeV', write_to_file=True, show_plot=False)
# calc_prod_rates_in_foil(isotope_list_list_Fe02, isotope_chain_parent_list_Fe02, 'Fe02', '50MeV', write_to_file=True, show_plot=False)
# calc_prod_rates_in_foil(isotope_list_list_Fe03, isotope_chain_parent_list_Fe03, 'Fe03', '50MeV', write_to_file=True, show_plot=False)
# calc_prod_rates_in_foil(isotope_list_list_Fe04, isotope_chain_parent_list_Fe04, 'Fe04', '50MeV', write_to_file=True, show_plot=False)
# calc_prod_rates_in_foil(isotope_list_list_Fe05, isotope_chain_parent_list_Fe05, 'Fe05', '50MeV', write_to_file=True, show_plot=False)



