
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import datetime


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



time_array = np.linspace(0,20,1000)*24*3600 #[s]
half_life_48V = 15.9735*24*3600 #[s], 48V
decay_constant = decay_constant_func(half_life_48V)



#Trying to impport real data, and not hard code dummy data:)
path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/CJ010317_Ti01_18cm_30MeV/'
file = 'CJ010317_Ti01_18cm_30MeV_peak_data.csv'

df = pd.read_csv(path+file,
        header=0,
        usecols=['isotope', 'decays', 'unc_decays', 'start_time', 'real_time'])

df['start_time'] = pd.to_datetime(df['start_time'])
df_48V = df[df['isotope'] == '48V']


# datetime(year, month, day, hour, minute, second)
time_EOB = datetime.datetime(2017, 2, 13, 12, 50, 3) #OBS: This is not the real end-of-beam time

test_delta_t = df_48V.iloc[0,3] - time_EOB
# print(test_delta_t.total_seconds())

decays = np.array(df_48V['decays'][:-1])
unc_decays = np.array(df_48V['unc_decays'][:-1])

start_times = []
for i in range(len(df_48V.iloc[:-1,3])):
    start_time_not_in_sec = df_48V.iloc[i,3] - time_EOB
    start_times.append(start_time_not_in_sec.total_seconds())

start_times = np.array(start_times)
real_time = np.array(df_48V.iloc[:-1,4])



#Calculating the A with uncertainty given the data I have. 
A_data, A_unc = calculate_A(decays, unc_decays, start_times, real_time)

print('A_data: ',A_data)
print('A_unc: ',A_unc)


#Fitting the activity curve to the data:
initial_A0_guess = 100
params, covariance = curve_fit(fit_function, start_times, A_data, sigma=A_unc, p0=[initial_A0_guess])
optimized_A0 = params[0]
activity_fit = activity_func(decay_constant, optimized_A0, time_array)


# Plot the original data and the fitted decay curve
plt.errorbar(start_times/(24*3600), A_data, yerr=A_unc, label='Data', color='lightgreen', fmt='o')
plt.plot(time_array/(24*3600), activity_fit, label=f'Fitted Decay Curve: N0={optimized_A0:.2f}, λ={decay_constant:.2f}', color='hotpink')
plt.xlabel('Time (d)')

plt.ylabel('Activity (Bq)')
plt.legend()
plt.show()




# # printing dataframe
# print('df:__________________________________ \n', df)
# print('df_48V:______________________________ \n', df_48V)




