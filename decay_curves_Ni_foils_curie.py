
import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import numpy as np 



# Ni peak data _____________________________________________________________________________________________________________________
path_Ni = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Ni_foils/'
file_AYNi02 = 'AY130217_Ni02_18cm_30MeV/AY130217_Ni02_18cm_30MeV_peak_data.csv'
file_AZNi02 = 'AZ130217_Ni02_18cm_30MeV/AZ130217_Ni02_18cm_30MeV_peak_data.csv'
file_BBNi01 = 'BB130217_Ni01_18cm_30MeV/BB130217_Ni01_18cm_30MeV_peak_data.csv'
file_BCNi03 = 'BC130217_Ni03_18cm_30MeV/BC130217_Ni03_18cm_30MeV_peak_data.csv'
file_BDNi04 = 'BD130217_Ni04_18cm_30MeV/BD130217_Ni04_18cm_30MeV_peak_data.csv'
file_BENi05 = 'BE130217_Ni05_18cm_30MeV/BE130217_Ni05_18cm_30MeV_peak_data.csv'
file_DANi01 = 'DA200317_Ni01_18cm_30MeV/DA200317_Ni01_18cm_30MeV_peak_data.csv'
file_DBNi02 = 'DB200317_Ni02_18cm_30MeV/DB200317_Ni02_18cm_30MeV_peak_data.csv'
file_DCNi03 = 'DC200317_Ni03_18cm_30MeV/DC200317_Ni03_18cm_30MeV_peak_data.csv'
file_DENi03 = 'DE210317_Ni03_18cm_30MeV/DE210317_Ni03_18cm_30MeV_peak_data.csv'
file_DFNi04 = 'DF240317_Ni04_18cm_30MeV/DF240317_Ni04_18cm_30MeV_peak_data.csv'
file_DGNi05 = 'DG240317_Ni05_18cm_30MeV/DG240317_Ni05_18cm_30MeV_peak_data.csv'


file_concat_Ni = 'combined_peak_data_Ni.csv'

df_AYNi02 = pd.read_csv(path_Ni+file_AYNi02)
df_AZNi02 = pd.read_csv(path_Ni+file_AZNi02)
df_BBNi01 = pd.read_csv(path_Ni+file_BBNi01)
df_BCNi03 = pd.read_csv(path_Ni+file_BCNi03)
df_BDNi04 = pd.read_csv(path_Ni+file_BDNi04)
df_BENi05 = pd.read_csv(path_Ni+file_BENi05)
df_DANi01 = pd.read_csv(path_Ni+file_DANi01)
df_DBNi02 = pd.read_csv(path_Ni+file_DBNi02)
df_DCNi03 = pd.read_csv(path_Ni+file_DCNi03)
df_DENi03 = pd.read_csv(path_Ni+file_DENi03)
df_DFNi04 = pd.read_csv(path_Ni+file_DFNi04)
df_DGNi05 = pd.read_csv(path_Ni+file_DGNi05)


# print(df_BENi05['isotope'])
# print(df_BENi05['decay_rate'])
# print(df_BDNi04['isotope'])

df_concat_Ni = pd.concat((df_AYNi02, df_AZNi02, df_BBNi01, df_BCNi03, df_BDNi04, df_BENi05, df_DANi01, df_DBNi02, df_DCNi03, df_DENi03, df_DFNi04, df_DGNi05), axis = 0)

# Specify the isotope and allowed energies
isotope_to_exclude = '58CO'
allowed_energies = [810.7593]

# Create a boolean mask based on the conditions
mask = (df_concat_Ni['isotope'] == isotope_to_exclude) & (~df_concat_Ni['energy'].isin(allowed_energies))

# Apply the mask to exclude rows
df_concat_Ni = df_concat_Ni[~mask]
df_concat_Ni = df_concat_Ni[~((df_concat_Ni['isotope'] == '58CO') & (df_concat_Ni['filename'].str.contains('30MeV/BB|30MeV/AY|30MeV/AZ|30MeV/BC|30MeV/BD|30MeV/BE')))]

# print(df_concat_Ni[df_concat_Ni['isotope']=='61CU'])

df_concat_Ni.to_csv(path_Ni+file_concat_Ni)



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
        isotopes, R, cov_R = fit_prod_rate(isotope_list, parent, foil, path_Ni, file_concat_Ni, stack, plot=show_plot)

        print(isotopes)
        print(R)
        # print(cov_R)
        if write_to_file==True:
            with open(csv_file_path, 'a', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                for i, isotope in enumerate(isotopes):
                    # print(i, isotope)
                    csv_writer.writerow([f'{isotope}', f'{R[i]}', f'{np.sqrt(cov_R[i,i])}', f'{np.sqrt(cov_R[i,i])/R[i]*100:.2e}'])




       

# isotope_list_list_Ni01 = [['52MN']]
# isotope_chain_parent_list_Ni01 = ['52MN']

# isotope_list_list_Ni02 = [['52MN']]
# isotope_chain_parent_list_Ni02 = ['52MN']



# isotope_list_list_Ni01 = [['54MN']]
# isotope_chain_parent_list_Ni01 = ['54MN']



# isotope_list_list_Ni01 = [['55CO']]
# isotope_chain_parent_list_Ni01 = ['55CO']

# isotope_list_list_Ni02 = [['55CO']]
# isotope_chain_parent_list_Ni02 = ['55CO']

# isotope_list_list_Ni03 = [['55CO']]
# isotope_chain_parent_list_Ni03 = ['55CO']

# isotope_list_list_Ni05 = [['55CO']]
# isotope_chain_parent_list_Ni05 = ['55CO']



# isotope_list_list_Ni01 = [['56CO']]
# isotope_chain_parent_list_Ni01 = ['56CO']

# isotope_list_list_Ni02 = [['56CO']]
# isotope_chain_parent_list_Ni02 = ['56CO']

# isotope_list_list_Ni03 = [['56CO']]
# isotope_chain_parent_list_Ni03 = ['56CO']

# isotope_list_list_Ni04 = [['56CO']]
# isotope_chain_parent_list_Ni04 = ['56CO']




# isotope_list_list_Ni01 = [['57NI', '57CO']]
# isotope_chain_parent_list_Ni01 = ['57NI']

# isotope_list_list_Ni02 = [['57NI', '57CO']]
# isotope_chain_parent_list_Ni02 = ['57NI']

# isotope_list_list_Ni03 = [['57NI', '57CO']]
# isotope_chain_parent_list_Ni03 = ['57NI']

# isotope_list_list_Ni04 = [['57NI', '57CO']]
# isotope_chain_parent_list_Ni04 = ['57NI']

# isotope_list_list_Ni05 = [['57NI', '57CO']]
# isotope_chain_parent_list_Ni05 = ['57NI']



# isotope_list_list_Ni01 = [['58CO']]
# isotope_chain_parent_list_Ni01 = ['58CO']

# isotope_list_list_Ni02 = [['58CO']]
# isotope_chain_parent_list_Ni02 = ['58CO']

# isotope_list_list_Ni03 = [['58CO']]
# isotope_chain_parent_list_Ni03 = ['58CO']

# isotope_list_list_Ni04 = [['58CO']]
# isotope_chain_parent_list_Ni04 = ['58CO']

# isotope_list_list_Ni05 = [['58CO']]
# isotope_chain_parent_list_Ni05 = ['58CO']



# isotope_list_list_Ni01 = [['60CO']]
# isotope_chain_parent_list_Ni01 = ['60CO']

# isotope_list_list_Ni04 = [['60CO']]
# isotope_chain_parent_list_Ni04 = ['60CO']

# isotope_list_list_Ni05 = [['60CO']]




# isotope_list_list_Ni02 = [['60CU']]
# isotope_chain_parent_list_Ni02 = ['60CU']

# isotope_list_list_Ni03 = [['60CU']]
# isotope_chain_parent_list_Ni03 = ['60CU']



isotope_list_list_Ni01 = [['61CU']]
isotope_chain_parent_list_Ni01 = ['61CU']

isotope_list_list_Ni02 = [['61CU']]
isotope_chain_parent_list_Ni02 = ['61CU']

isotope_list_list_Ni03 = [['61CU']]
isotope_chain_parent_list_Ni03 = ['61CU']

isotope_list_list_Ni04 = [['61CU']]
isotope_chain_parent_list_Ni04 = ['61CU']

isotope_list_list_Ni05 = [['61CU']]
isotope_chain_parent_list_Ni05 = ['61CU']



# isotope_list_list_Ni02 = [['64CU']]
# isotope_chain_parent_list_Ni02 = ['64CU']

# isotope_list_list_Ni03 = [['64CU']]
# isotope_chain_parent_list_Ni03 = ['64CU']

# isotope_list_list_Ni04 = [['64CU']]
# isotope_chain_parent_list_Ni04 = ['64CU']




# isotope_list_list_Ni01 = [['65NI']]
# isotope_chain_parent_list_Ni01 = ['65NI']

# isotope_list_list_Ni02 = [['65NI']]
# isotope_chain_parent_list_Ni02 = ['65NI']

# isotope_list_list_Ni03 = [['65NI']]
# isotope_chain_parent_list_Ni03 = ['65NI']

# isotope_list_list_Ni04 = [['65NI']]
# isotope_chain_parent_list_Ni04 = ['65NI']





calc_prod_rates_in_foil(isotope_list_list_Ni01, isotope_chain_parent_list_Ni01, 'Ni01', '30MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Ni02, isotope_chain_parent_list_Ni02, 'Ni02', '30MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Ni03, isotope_chain_parent_list_Ni03, 'Ni03', '30MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Ni04, isotope_chain_parent_list_Ni04, 'Ni04', '30MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Ni05, isotope_chain_parent_list_Ni05, 'Ni05', '30MeV', write_to_file=True, show_plot=True)





