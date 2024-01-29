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


def bateman_unstable_daughter(decay_const_p, decay_const_d,A0_p,A0_d,time):
    activity_d = decay_const_d*(A0_p*decay_const_p*(np.exp(-decay_const_p)+np.exp(-decay_const_d))/(decay_const_p-decay_const_d) + A0_d*np.exp(-decay_const_d*time))

    return activity_d


def fit_function(t, A0_p, A0_d):
#This function will be fitted to the data points. Only want to fit the A0s, and not decay constants
    return bateman_unstable_daughter(decay_constant_58COm, decay_constant_58CO, A0_p, A0_d, t)



time_array = np.linspace(0,40,1000)*24*3600 #[s]
half_life_58CO = 70.86*24*3600 #[s], 58CO
half_life_58COm = 9.10*3600 #[s], 58COm
decay_constant_58COm = decay_constant_func(half_life_58COm)
decay_constant_58CO = decay_constant_func(half_life_58CO)






def decay_chain(foil, isotope, decay_const_p, decay_const_d, path, file_concat):
    df = pd.read_csv(path+file_concat)

    #Sort out the rows containing the right foil and isotope
    mask = (df['filename'].str.contains(foil) & (df['isotope'] == isotope))
    df = df[mask]

    #Get the decays w unc 
    decays = np.array(df['decays'])
    unc_decays = np.array(df['unc_decays'])

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
        A_data_58CO, A_unc_58CO = calculate_A(decays, unc_decays, start_times_array, real_times)
        #Fitting the activity curve to the data:
        initial_A0_p_guess = 100
        initial_A0_d_guess = 100

        print('len(start_times_array): ', len(start_times_array))
        print('len(A_data_58CO): ', len(A_data_58CO))



        params, covariance = curve_fit(fit_function, start_times_array, A_data_58CO, sigma=A_unc_58CO, p0=[initial_A0_p_guess, initial_A0_d_guess])
        optimized_A0_p = params[0]
        optimized_A0_d = params[1]
        optimized_A0_unc = np.sqrt(covariance[1,1])
        activity_fit = bateman_unstable_daughter(decay_const_p, decay_const_d,optimized_A0_p,optimized_A0_d,time_array)

        # Plot the original data and the fitted decay curve
        plt.errorbar(start_times_array/(24*3600), A_data_58CO, yerr=A_unc_58CO, label='Data', color='skyblue', fmt='o', capsize = 5.0)
        plt.plot(time_array/(24*3600), activity_fit, label=f'Fitted Decay Curve for 58CO: A0={optimized_A0_d:.2f}', color='hotpink')
        plt.plot(time_array/(24*3600), activity_func(decay_constant_58COm, optimized_A0_p, time_array),label=f'Fitted Decay Curve for 58COm: A0={optimized_A0_p:.2f}', color='orchid')
        plt.xlabel('Time (d)')
        plt.ylabel('Activity (Bq)')
        plt.title(f'{foil}, {isotope}')
        plt.ylim([0,5000])
        plt.legend()
        plt.show()

        return optimized_A0_d, optimized_A0_unc




# Ni peak data _____________________________________________________________________________________________________________________
path_Ni = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Ni_foils/'
file_AYNi02 = 'AY130217_Ni02_18cm_30MeV/AY130217_Ni02_18cm_30MeV_peak_data.csv'
file_AZNi02 = 'AZ130217_Ni02_18cm_30MeV/AZ130217_Ni02_18cm_30MeV_peak_data.csv'
file_BBNi01 = 'BB130217_Ni01_18cm_30MeV/BB130217_Ni01_18cm_30MeV_peak_data.csv'
file_BCNi03 = 'BC130217_Ni03_18cm_30MeV/BC130217_Ni03_18cm_30MeV_peak_data.csv'
file_BDNi04 = 'BD130217_Ni04_18cm_30MeV/BD130217_Ni04_18cm_30MeV_peak_data.csv'
file_BENi05 = 'BE130217_Ni05_18cm_30MeV/BE130217_Ni05_18cm_30MeV_peak_data_by_hand.csv'
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
df_BENi05 = pd.read_csv(path_Ni+file_BENi05, delimiter=';')
df_DANi01 = pd.read_csv(path_Ni+file_DANi01)
df_DBNi02 = pd.read_csv(path_Ni+file_DBNi02)
df_DCNi03 = pd.read_csv(path_Ni+file_DCNi03)
df_DENi03 = pd.read_csv(path_Ni+file_DENi03)
df_DFNi04 = pd.read_csv(path_Ni+file_DFNi04)
df_DGNi05 = pd.read_csv(path_Ni+file_DGNi05)

df_concat_Ni = pd.concat((df_AYNi02, df_AZNi02, df_BBNi01, df_BCNi03, df_BDNi04, df_BENi05, df_DANi01, df_DBNi02, df_DCNi03, df_DENi03, df_DFNi04, df_DGNi05), axis = 0)


# Specify the isotope and allowed energies
isotope_to_exclude = '58CO'
# allowed_energies = [1037.843, 1238.288, 1771.357]
allowed_energies = [810.7593]


# Create a boolean mask based on the conditions
mask = (df_concat_Ni['isotope'] == isotope_to_exclude) & (~df_concat_Ni['energy'].isin(allowed_energies))

# Apply the mask to exclude rows
df_concat_Ni = df_concat_Ni[~mask]


df_concat_Ni.to_csv(path_Ni+file_concat_Ni)





#Running the code___________________________________________________________________
foil_list = ['Ni01', 'Ni02', 'Ni03', 'Ni04', 'Ni05']
isotope_list = ['58CO']


for foil in foil_list:
    print(f'__________FOIL {foil}___________')
    csv_file_path = f'./Calculated_A0/{foil}_A0_by_hand.csv'
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Isotope', 'A0', 'A0_unc'])
        for isotope in isotope_list:
            A0_d, A0_unc = decay_chain(foil, isotope, decay_constant_58COm, decay_constant_58CO, path_Ni, file_concat_Ni)
            # csv_writer.writerow([f'{isotope}', f'{A0}', f'{A0_unc}'])