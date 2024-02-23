

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from datetime import datetime
import csv


def decay_constant_func(half_life):
    #Calculate the decay constant given the half life of a nucleus
    return np.log(2)/half_life

def calculate_A(decays, unc_decays, start_times, real_time):
    #Calculate the activies from my data
    delta_t = real_time
    A = decays/delta_t
    A_unc = unc_decays/delta_t
    return A, A_unc

def activity_func(decay_constant, A0, t):
    #Returns the activity as a function of time
    return A0*np.exp(-decay_constant*t)


def fit_function(t, A0):
    #This function will be fitted to the data points. Only want to fit the A0, and not decay constants
    return activity_func(decay_constant, A0, t)


def decay_chain(foil, isotope, decay_constant, path, file_concat):
    df = pd.read_csv(path+file_concat)

    #Sort out the rows containing the right foil and isotope
    mask = (df['filename'].str.contains(foil) & (df['isotope'] == isotope))
    df = df[mask]
    # print(foil, isotope)
    # print('Before relative unc check \n')
    # print(df)

    # df = df[df['unc_decays'] / df['decays'] <= 1]
    # print('After relative unc check \n')
    # print(df)
    #Get the decays w unc 
    decays = np.array(df['decays'])
    unc_decays = np.array(df['unc_decays'])

    # print('relative uncertainty decays: ', unc_decays/decays)


    start_times_list = []
    # datetime(year, month, day, hour, minute, second)
    time_EOB = datetime(2017, 2, 13, 14, 27, 00) 
    start_times = df['start_time']
    datetime_start_times = [datetime.strptime(date_str, '%m/%d/%Y %H:%M:%S') for date_str in start_times]
    for start_time in datetime_start_times:
        start_time_in_sec = (start_time - time_EOB).total_seconds()
        start_times_list.append(start_time_in_sec)


    start_times_array = np.array(start_times_list)
    real_times = np.array(df['real_time'])

    bool = df[df['filename'].str.contains(foil) & (df['isotope'] == isotope)].shape[0] > 0
    if(bool):
        #Calculating the A with uncertainty given the data I have. 
        A_data, A_unc = calculate_A(decays, unc_decays, start_times_array, real_times)
        #Fitting the activity curve to the data:
        initial_A0_guess = 100
        params, covariance = curve_fit(fit_function, start_times_array, A_data, sigma=A_unc, p0=[initial_A0_guess])
        optimized_A0 = params[0]
        optimized_A0_unc = np.sqrt(covariance[0,0])
        activity_fit = activity_func(decay_constant, optimized_A0, time_array)

        # Plot the original data and the fitted decay curve
        plt.errorbar(start_times_array/(24*3600), A_data, yerr=A_unc, label='Data', color='skyblue', fmt='o', capsize = 5.0)
        plt.plot(time_array/(24*3600), activity_fit, label=f'Fitted Decay Curve: A0={optimized_A0:.2f}', color='hotpink')
        plt.xlabel('Time (d)')
        plt.ylabel('Activity (Bq)')
        plt.title(f'{foil}, {isotope}')
        plt.ylim(0)
        plt.legend()
        plt.show()

        return optimized_A0, optimized_A0_unc





time_array = np.linspace(0,40,1000)*24*3600 #[s]
half_life_56Co = 77.236*24*3600 #[s], 56Co
half_life_58Co = 70.86*24*3600 #[s], 58Co
half_life_90Nb = 14.60*3600 #[s]
decay_constant_56Co = decay_constant_func(half_life_56Co)
decay_constant_58Co = decay_constant_func(half_life_58Co)
decay_constant_90Nb = decay_constant_func(half_life_90Nb)




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


# Specify the isotope and allowed energies
isotope_to_exclude = '56CO'
# allowed_energies = [1037.843, 1238.288, 1771.357]
allowed_energies = [1037.843, 1238.288]


# Create a boolean mask based on the conditions
mask = (df_concat_Ni['isotope'] == isotope_to_exclude) & (~df_concat_Ni['energy'].isin(allowed_energies))

# Apply the mask to exclude rows
df_concat_Ni = df_concat_Ni[~mask]
df_concat_Ni.to_csv(path_Ni+file_concat_Ni)








# Zr peak data _____________________________________________________________________________________________________________________
path_Zr = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Zr_foils/'
file_AXZR05 = 'AX130217_Zr05_18cm_30MeV/AX130217_Zr05_18cm_30MeV_peak_data.csv'
file_BAZR01 = 'BA130217_Zr01_18cm_30MeV/BA130217_Zr01_18cm_30MeV_peak_data.csv'
file_BIZR04 = 'BI140217_Zr04_18cm_30MeV/BI140217_Zr04_18cm_30MeV_peak_data.csv'
file_BPZR01 = 'BP150217_Zr01_18cm_30MeV/BP150217_Zr01_18cm_30MeV_peak_data.csv'
file_BRZR02 = 'BR150217_Zr02_18cm_30MeV/BR150217_Zr02_18cm_30MeV_peak_data.csv'
file_BSZR05 = 'BS160217_Zr05_18cm_30MeV/BS160217_Zr05_18cm_30MeV_peak_data.csv'
file_BTZR03 = 'CH260217_Zr03_18cm_30MeV/CH260217_Zr03_18cm_30MeV_peak_data.csv'
file_BXZR01 = 'BX190217_Zr01_18cm_30MeV/BX190217_Zr01_18cm_30MeV_peak_data.csv'
file_CHZR03 = 'BT160217_Zr03_18cm_30MeV/BT160217_Zr03_18cm_30MeV_peak_data.csv'
file_CRZR04 = 'CR060317_Zr04_18cm_30MeV/CR060317_Zr04_18cm_30MeV_peak_data.csv'

file_BFZR01 = 'BF130217_Zr01_40cm_30MeV/BF130217_Zr01_40cm_30MeV_peak_data.csv'
file_BGZR02 = 'BG130217_Zr02_40cm_30MeV/BG130217_Zr02_40cm_30MeV_peak_data.csv'
file_BHZR03 = 'BH130217_Zr03_40cm_30MeV/BH130217_Zr03_40cm_30MeV_peak_data.csv'
file_BOZR04 = 'BO150217_Zr04_40cm_30MeV/BO150217_Zr04_40cm_30MeV_peak_data.csv'

file_CYZR02 = 'CY090317_Zr02_10cm_30MeV/CY090317_Zr02_10cm_30MeV_peak_data.csv'
file_CZZR05 = 'CZ120317_Zr05_10cm_30MeV/CZ120317_Zr05_10cm_30MeV_peak_data.csv'

file_concat_Zr = 'combined_peak_data_Zr.csv'

df_AXXR05 = pd.read_csv(path_Zr+file_AXZR05)
df_BAZR01 = pd.read_csv(path_Zr+file_BAZR01)
df_BIZR04 = pd.read_csv(path_Zr+file_BIZR04)
df_BPZR01 = pd.read_csv(path_Zr+file_BPZR01)
df_BRZR02 = pd.read_csv(path_Zr+file_BRZR02)
df_BSZR05 = pd.read_csv(path_Zr+file_BSZR05)
df_BTZR03 = pd.read_csv(path_Zr+file_BTZR03)
df_BXZR01 = pd.read_csv(path_Zr+file_BXZR01)
df_CHZR03 = pd.read_csv(path_Zr+file_CHZR03)
df_CRZR04 = pd.read_csv(path_Zr+file_CRZR04)

df_BFZR01 = pd.read_csv(path_Zr+file_BFZR01)
df_BGZR02 = pd.read_csv(path_Zr+file_BGZR02)
df_BHZR03 = pd.read_csv(path_Zr+file_BHZR03)
df_BOZR04 = pd.read_csv(path_Zr+file_BOZR04)

df_CYZR02 = pd.read_csv(path_Zr+file_CYZR02)
df_CZZR05 = pd.read_csv(path_Zr+file_CZZR05)



df_concat_Zr = pd.concat((df_AXXR05, df_BAZR01, df_BIZR04, df_BPZR01, df_BRZR02, df_BSZR05, df_BTZR03, df_BXZR01, df_CHZR03, df_CRZR04, df_BFZR01, df_BGZR02, df_BHZR03, df_BOZR04, df_CYZR02, df_CZZR05), axis = 0)
df_concat_Zr= df_concat_Zr[(df_concat_Zr['isotope'] != '90NB') | (df_concat_Zr['energy'] != 329.058)] # Exclude rows where isotope is '90NB' and energy is 329.058



# isotopes_Zr = ['96NB', '90NB']
# isotopes_Zr = ['90NB']

# isotopes_Zr = ['87NB','89NB', '90NB', '95NB']
# isotopes_Zr = ['87NB']

# isotopes_Zr = ['88Y', '91Y', '92Y']


# foil_list_Zr = ['Zr01', 'Zr02', 'Zr03', 'Zr04', 'Zr05']
# foil_list_Zr = ['Zr04']







#Running the code______________________________________________________________________________________
# decay_chain(decay_constant, df_48V)
foil_list_Ni = ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05']
# isotope_list = ['56CO', '58CO', '61CU']
isotopes_Ni = ['56CO']

foil_list = foil_list_Ni
isotope_list = isotopes_Ni
path = path_Ni
file_concat = file_concat_Ni

for foil in foil_list:
    csv_file_path = f'./Calculated_A0/{foil}_A0_by_hand.csv'
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Isotope', 'A0', 'A0_unc'])
        for isotope in isotope_list:
            if isotope == '56CO':
                decay_constant = decay_constant_56Co
            elif isotope == '58CO':
                decay_constant = decay_constant_58Co
            elif isotope == '90NB':
                decay_constant = decay_constant_90Nb
            A0, A0_unc = decay_chain(foil, isotope, decay_constant, path, file_concat)
            csv_writer.writerow([f'{isotope}', f'{A0}', f'{A0_unc}'])