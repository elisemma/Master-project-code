
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def decay_constant_func(half_life):
    return np.log(2)/half_life

def calculate_N0(decay_constant, decays, unc_decays, start_times, stop_times):
    delta_t = stop_times-start_times
    A0 = decays/delta_t
    A0_unc = unc_decays/delta_t
    N0 = A0/decay_constant
    N0_unc = A0_unc/decay_constant
    return N0, N0_unc

def decay_curve_func(t, N0, decay_constant):
    return N0 * np.exp(-decay_constant * t)

def fit_function(t, N0):
    return decay_curve_func(t, N0, decay_constant)



time_array = np.linspace(0,8,1000)
half_life = 15.9735 #[d], 48V
decay_constant = decay_constant_func(half_life)

decays = np.array([6.0e5, 7.0e5])
unc_decays = np.array([2.0e4, 3.0e4])
start_times = np.array([5.0, 6.0])
stop_times = np.array([5.1, 6.1])

N0_data, N0_unc = calculate_N0(decay_constant, decays, unc_decays, start_times,stop_times)

initial_N0_guess = 1000
params, covariance = curve_fit(fit_function, start_times, N0_data, sigma=N0_unc, p0=[initial_N0_guess])
optimized_N0 = params[0]
activity_fit = decay_curve_func(time_array, optimized_N0, decay_constant)



# Plot the original data and the fitted decay curve
plt.errorbar(start_times, N0_data, yerr=N0_unc, label='Data', color='hotpink', fmt='o')
plt.plot(time_array, activity_fit, label=f'Fitted Decay Curve: N0={optimized_N0:.2f}, λ={decay_constant:.2f}', color='hotpink')
plt.xlabel('Time (d)')
plt.ylabel('Activity (Bq)')
plt.legend()
plt.show()



