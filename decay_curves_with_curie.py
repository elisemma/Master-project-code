import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def fit_prod_rate(isotope, foil, path, file):
    R_estimated = 5
    t_irr_h = 0.33
    dc = ci.DecayChain(isotope, units='h', R=[[R_estimated, t_irr_h]])
    dc.get_counts(foil, '02/13/2017 14:27:00', path+file)

    isotopes, r, var_r = dc.fit_R()
    R = float(r[0])    
    var_R = float(var_r[0][0])
    dc.plot(palette='american', shade='dark', style='paper')
    return R, var_R



def generate_prod_rate_csv(foil_type, isotopes, foil_list, path, file_concat):

    R_file = open(path+'production_rate_'+foil_type+'.csv', 'w')
    R_file.write('foil, isotope, R, var_R \n')
    print('foil, isotope, R, var_R \n')

    df = pd.read_csv(path+file_concat)
    for isotope in isotopes:
        for foil in foil_list:
            bool = df[df['filename'].str.contains(foil) & (df['isotope'] == isotope)].shape[0] > 0
            if(bool):
                print(f'Foil: {foil}_________________________')
                R, var_R = fit_prod_rate(isotope, foil, path, file_concat)
                print(f'{foil}, {isotope}, {R}, {var_R} \n')
                R_file.write(f'{foil}, {isotope}, {R}, {var_R} \n')
    R_file.close()




# Ti peak data _____________________________________________________________________________________________________________________
path_Ti = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Ti_foils/'
file_CJTi01 = 'CJ010317_Ti01_18cm_30MeV/CJ010317_Ti01_18cm_30MeV_peak_data.csv'
file_CKTi02 = 'CK010317_Ti02_18cm_30MeV/CK010317_Ti02_18cm_30MeV_peak_data.csv'
file_CLTi03 = 'CL010317_Ti03_18cm_30MeV/CL010317_Ti03_18cm_30MeV_peak_data.csv'
file_CMTi04 = 'CM010317_Ti04_18cm_30MeV/CM010317_Ti04_18cm_30MeV_peak_data.csv'
file_CPTi05 = 'CP030317_Ti05_18cm_30MeV/CP030317_Ti05_18cm_30MeV_peak_data.csv'
file_CQTi04 = 'CQ030317_Ti04_18cm_30MeV/CQ030317_Ti04_18cm_30MeV_peak_data.csv'
file_CSTi01 = 'CS060317_Ti01_18cm_30MeV/CS060317_Ti01_18cm_30MeV_peak_data.csv'
file_CTTi02 = 'CT060317_Ti02_18cm_30MeV/CT060317_Ti02_18cm_30MeV_peak_data.csv'
file_CUTi03 = 'CU060317_Ti03_18cm_30MeV/CU060317_Ti03_18cm_30MeV_peak_data.csv'

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


#Want to exclude one bad peak:
# Values to exclude
isotope_to_exclude = '48V'
energy_to_exclude = 928.327
# Create a boolean mask based on the conditions
mask = (df_CKTi02['isotope'] != isotope_to_exclude) | (df_CKTi02['energy'] != energy_to_exclude)
# Use boolean indexing to exclude rows based on the conditions
df_CKTi02 = df_CKTi02[mask]


df_concat_Ti = pd.concat((df_CJTi01, df_CKTi02, df_CLTi03, df_CMTi04, df_CPTi05, df_CQTi04, df_CSTi01, df_CTTi02, df_CUTi03), axis = 0)
df_concat_Ti.to_csv(path_Ti+file_concat_Ti)

# isotopes_Ti = ['48V', '46SC']
# isotopes_Ti = ['48V']

# foil_list_Ti = ['Ti01', 'Ti02','Ti03', 'Ti04', 'Ti05']
# foil_list_Ti = ['Ti02']




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


df_concat_Ni = pd.concat((df_AYNi02, df_AZNi02, df_BBNi01, df_BCNi03, df_BDNi04, df_BENi05, df_DANi01, df_DBNi02, df_DCNi03, df_DENi03, df_DFNi04, df_DGNi05), axis = 0)
df_concat_Ni.to_csv(path_Ni+file_concat_Ni)

# isotopes_Ni = ['56CO', '58CO', '61CU']
isotopes_Ni = ['56CO']

# foil_list_Ni = ['Ni01', 'Ni02','Ni03', 'Ni04', 'Ni05']
foil_list_Ni = ['Ni01']




#Finding huge Chi2:
# print(df_concat_Ti[df_concat_Ti['chi2']>10])
# print(df_concat_Ni[df_concat_Ni['chi2']>10])

# print(df_concat_Ti[df_concat_Ti['filename'].str.contains('Ti01') & (df_concat_Ti['isotope'] == '48V')])
df_concat_Ni = df_concat_Ni.drop(columns=['start_time'])
df_concat_Ni = df_concat_Ni.drop(columns=['live_time'])
df_concat_Ni = df_concat_Ni.drop(columns=['real_time'])
df_concat_Ni = df_concat_Ni.drop(columns=['intensity'])
df_concat_Ni = df_concat_Ni.drop(columns=['unc_efficiency'])
df_concat_Ni = df_concat_Ni.drop(columns=['efficiency'])





# print(df_concat_Ti[df_concat_Ti['filename'].str.contains('Ti02') & (df_concat_Ti['isotope'] == '48V')])
print(df_concat_Ni[df_concat_Ni['filename'].str.contains('Ni01') & (df_concat_Ni['isotope'] == '56CO')])



#Running the code_______________________________________________________________________________________________________________________
# generate_prod_rate_csv('Ti', isotopes_Ti, foil_list_Ti, path_Ti, file_concat_Ti)
generate_prod_rate_csv('Ni', isotopes_Ni, foil_list_Ni, path_Ni, file_concat_Ni)




