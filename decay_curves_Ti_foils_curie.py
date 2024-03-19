
import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import numpy as np 



# Ti peak data _____________________________________________________________________________________________________________________
path_Ti = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Ti_foils/'
#30 MeV
file_CJTi01 = 'CJ010317_Ti01_18cm_30MeV/CJ010317_Ti01_18cm_30MeV_peak_data.csv'
file_CKTi02 = 'CK010317_Ti02_18cm_30MeV/CK010317_Ti02_18cm_30MeV_peak_data.csv'
file_CLTi03 = 'CL010317_Ti03_18cm_30MeV/CL010317_Ti03_18cm_30MeV_peak_data.csv'
file_CMTi04 = 'CM010317_Ti04_18cm_30MeV/CM010317_Ti04_18cm_30MeV_peak_data.csv'
file_CPTi05 = 'CP030317_Ti05_18cm_30MeV/CP030317_Ti05_18cm_30MeV_peak_data.csv'
file_CQTi04 = 'CQ030317_Ti04_18cm_30MeV/CQ030317_Ti04_18cm_30MeV_peak_data.csv'
file_CSTi01 = 'CS060317_Ti01_18cm_30MeV/CS060317_Ti01_18cm_30MeV_peak_data.csv'
file_CTTi02 = 'CT060317_Ti02_18cm_30MeV/CT060317_Ti02_18cm_30MeV_peak_data.csv'
file_CUTi03 = 'CU060317_Ti03_18cm_30MeV/CU060317_Ti03_18cm_30MeV_peak_data.csv'
#50 MeV
file_CCTi06 = 'CC220217_Ti06_18cm_50MeV/CC220217_Ti06_18cm_50MeV_peak_data.csv'
file_CETi08 = 'CE230217_Ti08_18cm_50MeV/CE230217_Ti08_18cm_50MeV_peak_data.csv'
file_CFTi09 = 'CF240217_Ti09_18cm_50MeV/CF240217_Ti09_18cm_50MeV_peak_data.csv'
file_CGTi10 = 'CG240217_Ti10_18cm_50MeV/CG240217_Ti10_18cm_50MeV_peak_data.csv'
file_CITi11 = 'CI010317_Ti11_18cm_50MeV/CI010317_Ti11_18cm_50MeV_peak_data.csv'
file_COTi06 = 'CO020317_Ti06_18cm_50MeV/CO020317_Ti06_18cm_50MeV_peak_data.csv'
# file_CWTi11 = 'CW080317_Ti11_10cm_50MeV/CW080317_Ti11_10cm_50MeV_peak_data.csv'

file_concat_Ti = 'combined_peak_data_Ti.csv'

df_CJTi01 = pd.read_csv(path_Ti+file_CJTi01)
df_CKTi02 = pd.read_csv(path_Ti+file_CKTi02)
df_CLTi03 = pd.read_csv(path_Ti+file_CLTi03)
df_CMTi04 = pd.read_csv(path_Ti+file_CMTi04)
df_CPTi05 = pd.read_csv(path_Ti+file_CPTi05)
df_CQTi04 = pd.read_csv(path_Ti+file_CQTi04)
df_CSTi01 = pd.read_csv(path_Ti+file_CSTi01)
df_CTTi02 = pd.read_csv(path_Ti+file_CTTi02)
df_CUTi03 = pd.read_csv(path_Ti+file_CUTi03)

df_CCTi06 = pd.read_csv(path_Ti+file_CCTi06)
df_CETi08 = pd.read_csv(path_Ti+file_CETi08)
df_CFTi09 = pd.read_csv(path_Ti+file_CFTi09)
df_CGTi10 = pd.read_csv(path_Ti+file_CGTi10)
df_CITi11 = pd.read_csv(path_Ti+file_CITi11)
df_COTi06 = pd.read_csv(path_Ti+file_COTi06)
# df_CWTi11 = pd.read_csv(path_Ti+file_CWTi11)



df_concat_Ti = pd.concat((df_CJTi01, df_CKTi02, df_CLTi03, df_CMTi04, df_CPTi05, df_CQTi04, df_CSTi01, df_CTTi02, df_CUTi03, df_CCTi06, df_CETi08, df_CFTi09, df_CGTi10, df_CITi11, df_COTi06), axis = 0)
df_concat_Ti = df_concat_Ti[~((df_concat_Ti['isotope'] == '48V') & (~df_concat_Ti['energy'].isin([944.130, 983.525, 1312.106])))]

df_concat_Ti.to_csv(path_Ti+file_concat_Ti)





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
        isotopes, R, cov_R = fit_prod_rate(isotope_list, parent, foil, path_Ti, file_concat_Ti, stack, plot=show_plot)

        print(isotopes)
        print(R)
        # print(cov_R)
        if write_to_file==True:
            with open(csv_file_path, 'a', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                for i, isotope in enumerate(isotopes):
                    # print(i, isotope)
                    csv_writer.writerow([f'{isotope}', f'{R[i]}', f'{np.sqrt(cov_R[i,i])}', f'{np.sqrt(cov_R[i,i])/R[i]*100:.2e}'])




       


# isotope_list_list_Ti06 = [['44SC']]
# isotope_chain_parent_list_Ti06 = ['44SC']

# isotope_list_list_Ti08 = [['44SC']]
# isotope_chain_parent_list_Ti08 = ['44SC']

# isotope_list_list_Ti09 = [['44SC']]
# isotope_chain_parent_list_Ti09 = ['44SC']

# isotope_list_list_Ti10 = [['44SC']]
# isotope_chain_parent_list_Ti10 = ['44SC']





# isotope_list_list_Ti01 = [['46SC']]
# isotope_chain_parent_list_Ti01 = ['46SC']

# isotope_list_list_Ti02 = [['46SC']]
# isotope_chain_parent_list_Ti02 = ['46SC']

# isotope_list_list_Ti03 = [['46SC']]
# isotope_chain_parent_list_Ti03 = ['46SC']

# isotope_list_list_Ti04 = [['46SC']]
# isotope_chain_parent_list_Ti04 = ['46SC']

# isotope_list_list_Ti05 = [['46SC']]
# isotope_chain_parent_list_Ti05 = ['46SC']

# isotope_list_list_Ti06 = [['46SC']]
# isotope_chain_parent_list_Ti06 = ['46SC']

# isotope_list_list_Ti08 = [['46SC']]
# isotope_chain_parent_list_Ti08 = ['46SC']

# isotope_list_list_Ti09 = [['46SC']]
# isotope_chain_parent_list_Ti09 = ['46SC']

# isotope_list_list_Ti10 = [['46SC']]
# isotope_chain_parent_list_Ti10 = ['46SC']

# isotope_list_list_Ti11 = [['46SC']]
# isotope_chain_parent_list_Ti11 = ['46SC']




# isotope_list_list_Ti01 = [['47SC']]
# isotope_chain_parent_list_Ti01 = ['47SC']

# isotope_list_list_Ti02 = [['47SC']]
# isotope_chain_parent_list_Ti02 = ['47SC']

# isotope_list_list_Ti03 = [['47SC']]
# isotope_chain_parent_list_Ti03 = ['47SC']

# isotope_list_list_Ti06 = [['47SC']]
# isotope_chain_parent_list_Ti06 = ['47SC']

# isotope_list_list_Ti08 = [['47SC']]
# isotope_chain_parent_list_Ti08 = ['47SC']

# isotope_list_list_Ti09 = [['47SC']]
# isotope_chain_parent_list_Ti09 = ['47SC']

# isotope_list_list_Ti10 = [['47SC']]
# isotope_chain_parent_list_Ti10 = ['47SC']

# isotope_list_list_Ti11 = [['47SC']]
# isotope_chain_parent_list_Ti11 = ['47SC']



# isotope_list_list_Ti01 = [['48SC']]
# isotope_chain_parent_list_Ti01 = ['48SC']

# isotope_list_list_Ti02 = [['48SC']]
# isotope_chain_parent_list_Ti02 = ['48SC']

# isotope_list_list_Ti03 = [['48SC']]
# isotope_chain_parent_list_Ti03 = ['48SC']

# isotope_list_list_Ti04 = [['48SC']]
# isotope_chain_parent_list_Ti04 = ['48SC']

# isotope_list_list_Ti06 = [['48SC']]
# isotope_chain_parent_list_Ti06 = ['48SC']

# isotope_list_list_Ti08 = [['48SC']]
# isotope_chain_parent_list_Ti08 = ['48SC']

# isotope_list_list_Ti09 = [['48SC']]
# isotope_chain_parent_list_Ti09 = ['48SC']

# isotope_list_list_Ti10 = [['48SC']]
# isotope_chain_parent_list_Ti10 = ['48SC']

# isotope_list_list_Ti11 = [['48SC']]
# isotope_chain_parent_list_Ti11 = ['48SC']



# isotope_list_list_Ti01 = [['48V']]
# isotope_chain_parent_list_Ti01 = ['48V']

# isotope_list_list_Ti02 = [['48V']]
# isotope_chain_parent_list_Ti02 = ['48V']

# isotope_list_list_Ti03 = [['48V']]
# isotope_chain_parent_list_Ti03 = ['48V']

# isotope_list_list_Ti06 = [['48V']]
# isotope_chain_parent_list_Ti06 = ['48V']

# isotope_list_list_Ti08 = [['48V']]
# isotope_chain_parent_list_Ti08 = ['48V']

# isotope_list_list_Ti09 = [['48V']]
# isotope_chain_parent_list_Ti09 = ['48V']

# isotope_list_list_Ti10 = [['48V']]
# isotope_chain_parent_list_Ti10 = ['48V']

# isotope_list_list_Ti11 = [['48V']]
# isotope_chain_parent_list_Ti11 = ['48V']




isotope_list_list_Ti01 = [['46SC'], ['47SC'], ['48V']]
isotope_chain_parent_list_Ti01 = ['46SC', '47SC', '48V']

isotope_list_list_Ti02 = [['46SC'], ['47SC'], ['48V']]
isotope_chain_parent_list_Ti02 = ['46SC', '47SC', '48V']

isotope_list_list_Ti03 = [['46SC'], ['47SC'], ['48V']]
isotope_chain_parent_list_Ti03 = ['46SC', '47SC', '48V']

isotope_list_list_Ti04 = [['46SC'], ['48V']]
isotope_chain_parent_list_Ti04 = ['46SC', '48V']

isotope_list_list_Ti05 = [['46SC'], ['48V']]
isotope_chain_parent_list_Ti05 = ['46SC', '48V']

isotope_list_list_Ti06 = [['44SCm', '44SC'], ['46SC'], ['47SC'],['48SC'], ['48V']]
isotope_chain_parent_list_Ti06 = ['44SCm', '46SC', '47SC', '48SC', '48V']

isotope_list_list_Ti08 = [['44SCm', '44SC'], ['46SC'], ['47SC'],['48SC'], ['48V']]
isotope_chain_parent_list_Ti08 = ['44SCm', '46SC', '47SC', '48SC', '48V']

isotope_list_list_Ti09 = [['44SCm', '44SC'], ['46SC'], ['47SC'],['48SC'], ['48V']]
isotope_chain_parent_list_Ti09 = ['44SCm', '46SC', '47SC', '48SC', '48V']

isotope_list_list_Ti10 = [['44SCm', '44SC'], ['46SC'], ['47SC'],['48SC'], ['48V']]
isotope_chain_parent_list_Ti10 = ['44SCm', '46SC', '47SC', '48SC', '48V']

isotope_list_list_Ti11 = [['46SC'], ['47SC'], ['48V']]
isotope_chain_parent_list_Ti11 = ['46SC', '47SC', '48V']










calc_prod_rates_in_foil(isotope_list_list_Ti01, isotope_chain_parent_list_Ti01, 'Ti01', '30MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti02, isotope_chain_parent_list_Ti02, 'Ti02', '30MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti03, isotope_chain_parent_list_Ti03, 'Ti03', '30MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti04, isotope_chain_parent_list_Ti04, 'Ti04', '30MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti05, isotope_chain_parent_list_Ti05, 'Ti05', '30MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti06, isotope_chain_parent_list_Ti06, 'Ti06', '50MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti08, isotope_chain_parent_list_Ti08, 'Ti08', '50MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti09, isotope_chain_parent_list_Ti09, 'Ti09', '50MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti10, isotope_chain_parent_list_Ti10, 'Ti10', '50MeV', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Ti11, isotope_chain_parent_list_Ti11, 'Ti11', '50MeV', write_to_file=True, show_plot=False)




