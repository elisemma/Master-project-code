import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

def fit_activity(isotope, foil, path, file):
    t_irr_h = 0.33
    dc = ci.DecayChain(isotope, units='h', A0=1000)
    dc.get_counts(foil, '02/13/2017 14:27:00', path+file)

    isotopes, a0, var_a0 = dc.fit_A0()
    
    if isotope == '58COm':
        print('_______________')
        print(isotopes)
        A0 = float(a0[1])
        var_A0 = float(var_a0[1][1])
    else:
        A0 = float(a0[0])
        var_A0 = float(var_a0[0][0])

    unc_A0 = np.sqrt(var_A0)
    dc.plot(palette='american', shade='dark', style='paper', max_plot_chi2=100)
    return A0, unc_A0



def generate_activity_csv(isotope_list, foil_list, path, file_concat):
    df = pd.read_csv(path+file_concat)
    for foil in foil_list:
        print(foil)
        print(df[df['filename'].str.contains(foil) & (df['isotope'] == '56CO')])

        csv_file_path = f'./Calculated_A0/{foil}_A0_by_curie.csv'
        with open(csv_file_path, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['Isotope', 'A0', 'A0_unc'])
            for isotope in isotope_list:
                bool = df[df['filename'].str.contains(foil) & (df['isotope'] == isotope)].shape[0] > 0
                if(bool):
                    A0, A0_unc = fit_activity(isotope, foil, path, file_concat)
                    csv_writer.writerow([f'{isotope}', f'{A0}', f'{A0_unc}'])
                   



# Ti peak data _____________________________________________________________________________________________________________________
path_Ti = './MyGeneratedFiles/Ti_foils/'
file_CCTi06 = 'CC220217_Ti06_18cm_50MeV/CC220217_Ti06_18cm_50MeV_peak_data.csv'
file_CETi08 = 'CE230217_Ti08_18cm_50MeV/CE230217_Ti08_18cm_50MeV_peak_data.csv'
file_CFTi09 = 'CF240217_Ti09_18cm_50MeV/CF240217_Ti09_18cm_50MeV_peak_data.csv'
file_CGTi10 = 'CG240217_Ti10_18cm_50MeV/CG240217_Ti10_18cm_50MeV_peak_data.csv'
file_CITi11 = 'CI010317_Ti11_18cm_50MeV/CI010317_Ti11_18cm_50MeV_peak_data.csv'
file_COTi06 = 'CO020317_Ti06_18cm_50MeV/CO020317_Ti06_18cm_50MeV_peak_data.csv'
file_CWTi11 = 'CW080317_Ti11_10cm_50MeV/CW080317_Ti11_10cm_50MeV_peak_data.csv'

file_concat_Ti = 'combined_peak_data_Ti.csv'


df_CCTi06 = pd.read_csv(path_Ti+file_CCTi06)
df_CETi08 = pd.read_csv(path_Ti+file_CETi08)
df_CFTi09 = pd.read_csv(path_Ti+file_CFTi09)
df_CGTi10 = pd.read_csv(path_Ti+file_CGTi10)
df_CITi11 = pd.read_csv(path_Ti+file_CITi11)
df_COTi06 = pd.read_csv(path_Ti+file_COTi06)
df_CWTi11 = pd.read_csv(path_Ti+file_CWTi11)


df_concat_Ti = pd.concat((df_CCTi06, df_CETi08, df_CFTi09, df_CGTi10, df_CITi11, df_COTi06, df_CWTi11), axis = 0)
df_concat_Ti = df_concat_Ti[~((df_concat_Ti['isotope'] == '48V') & (~df_concat_Ti['energy'].isin([944.130, 983.525, 1312.106])))]
df_concat_Ti.to_csv(path_Ti+file_concat_Ti)

isotopes_Ti = ['46SC', '48V']

foil_list_Ti = ['Ti06', 'Ti08','Ti09', 'Ti10', 'Ti11']


# print(df_concat_Ti[df_concat_Ti['filename'].str.contains('Ti06') & (df_concat_Ti['isotope'] == '46SC')])






# Fe peak data _____________________________________________________________________________________________________________________
path_Fe = './MyGeneratedFiles/Fe_foils/'
file_BUFe01 = 'BU170217_Fe01_18cm_50MeV/BU170217_Fe01_18cm_50MeV_peak_data.csv'
file_BVFe02 = 'BV170217_Fe02_18cm_50MeV/BV170217_Fe02_18cm_50MeV_peak_data.csv'
file_BWFe03 = 'BW170217_Fe03_18cm_50MeV/BW170217_Fe03_18cm_50MeV_peak_data.csv'
file_BYFe04 = 'BY190217_Fe04_18cm_50MeV/BY190217_Fe04_18cm_50MeV_peak_data.csv'
file_BZFe05 = 'BZ190217_Fe05_18cm_50MeV/BZ190217_Fe05_18cm_50MeV_peak_data.csv'

file_concat_Fe = 'combined_peak_data_Fe.csv'


df_BUFe01 = pd.read_csv(path_Fe+file_BUFe01)
df_BVFe02 = pd.read_csv(path_Fe+file_BVFe02)
df_BWFe03 = pd.read_csv(path_Fe+file_BWFe03)
df_BYFe04 = pd.read_csv(path_Fe+file_BYFe04)
df_BZFe05 = pd.read_csv(path_Fe+file_BZFe05)


df_concat_Fe = pd.concat((df_BUFe01, df_BVFe02, df_BWFe03, df_BYFe04, df_BZFe05), axis = 0)
df_concat_Fe = df_concat_Fe[~((df_concat_Fe['isotope'] == '56CO') & (~df_concat_Fe['energy'].isin([846.770, 1037.843, 1238.288, 1771.357, 2034.791, 2598.500])))]
# df_concat_Fe = df_concat_Fe[~((df_concat_Fe['isotope'] == '56CO') & (~df_concat_Fe['energy'].isin([1037.843, 1238.288, 1771.357, 2034.791, 2598.500])))]
df_concat_Fe.to_csv(path_Fe+file_concat_Fe)


isotopes_Fe = ['56CO']

foil_list_Fe = ['Fe01', 'Fe02','Fe03', 'Fe04', 'Fe05']











#Running the code_______________________________________________________________________________________________________________________
# generate_activity_csv(isotopes_Ti, foil_list_Ti, path_Ti, file_concat_Ti)
generate_activity_csv(isotopes_Fe, foil_list_Fe, path_Fe, file_concat_Fe)







