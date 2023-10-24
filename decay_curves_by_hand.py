
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def decay_constant_func(half_life):
    return np.log(2)/half_life

def calc_A0(decay_constant, decays, unc_decays, start_time, stop_time):
    count_time = stop_time - start_time
    after_EOB_time = start_time
    A0 = decays*decay_constant/(1-np.exp(-decay_constant*count_time))*np.exp(-decay_constant*after_EOB_time)
    A0_unc = unc_decays*decay_constant/(1-np.exp(-decay_constant*count_time))*np.exp(-decay_constant*after_EOB_time)
    return A0, A0_unc

def activity_func(decay_constant, A0, t):
    return A0*np.exp(-decay_constant*t)

def fit_function(t, A0):
    return activity_func(decay_constant, A0, t)



time_array = np.linspace(0,8,1000)*24*3600
half_life = 15.9735*24*3600 #[sec], 48V
decay_constant = decay_constant_func(half_life)

decays = np.array([6.0e5, 7.0e5])
unc_decays = np.array([2.0e4, 3.0e4])
start_times = np.array([5.0, 6.0])*24*3600
stop_times = np.array([5.1, 6.1])*24*3600

A0_data, A0_unc = calc_A0(decay_constant, decays, unc_decays, start_times,stop_times)

initial_A0_guess = 100
params, covariance = curve_fit(fit_function, start_times, A0_data, sigma=A0_unc, p0=[initial_A0_guess])
optimized_A0 = params[0]
activity_fit = activity_func(decay_constant, optimized_A0, time_array)



# Plot the original data and the fitted decay curve
plt.errorbar(start_times/(24*3600), A0_data, yerr=A0_unc, label='Data', color='hotpink', fmt='o')
plt.plot(time_array/(24*3600), activity_fit, label=f'Fitted Decay Curve: N0={optimized_A0:.2f}, λ={decay_constant:.2f}', color='hotpink')
plt.xlabel('Time (d)')
plt.ylabel('Activity (Bq)')
plt.legend()
plt.show()



