import numpy as np 
from datetime import datetime




def calc_decay_constant(t_half):
    return np.log(2)/t_half

def calc_delta_decay_const(t_half, delta_t_half):
    delta_decay_const = np.log(2)*delta_t_half/(t_half**2)
    return delta_decay_const

def calc_delta_A0(A0, A_t, delta_A_t, decay_const, delta_decay_const, t, delta_t):
    delta_A0 = A0*np.sqrt((delta_A_t/A_t)**2 + decay_const**2*t**2*((delta_decay_const/decay_const)**2+(delta_t/t)**2))
    return delta_A0


t_half_56Co = 77.236*24*3600 #[s]
delta_t_half_56Co = 0.026*24*3600 #[s]
delta_decay_const_56Co = calc_delta_decay_const(t_half_56Co, delta_t_half_56Co)
decay_const_56Co = calc_decay_constant(t_half_56Co)


time_EOB = datetime(2017, 2, 12, 19, 22, 00)
start_time_Fe01 = datetime(2017, 2, 17, 10, 44, 6) # 02/17/2017 10:44:06
start_time_Fe02 = datetime(2017, 2, 17, 13, 53, 4) # 02/17/2017 13:53:04
start_time_Fe03 = datetime(2017, 2, 17, 16, 21, 12) # 02/17/2017 16:21:12
start_time_Fe04 = datetime(2017, 2, 19, 11, 23, 38) # 02/19/2017 11:23:38
start_time_Fe05 = datetime(2017, 2, 19, 12, 52, 57) # 02/19/2017 12:52:57 



start_time_in_sec_Fe01 = (start_time_Fe01 - time_EOB).total_seconds()
start_time_in_sec_Fe02 = (start_time_Fe02 - time_EOB).total_seconds()
start_time_in_sec_Fe03 = (start_time_Fe03 - time_EOB).total_seconds()
start_time_in_sec_Fe04 = (start_time_Fe04 - time_EOB).total_seconds()
start_time_in_sec_Fe05 = (start_time_Fe05 - time_EOB).total_seconds()


#                      calc_delta_A0(A0,                 A_t,          delta_A_t, decay_const,       delta_decay_const,     t,                delta_t)
delta_A0_56Co_foil_1 = calc_delta_A0(655.9564121176588, 633.336787, 44.732989, decay_const_56Co, delta_decay_const_56Co, start_time_in_sec_Fe01, 2)
delta_A0_56Co_foil_2 = calc_delta_A0(861.0750687694533, 830.5102,  59.424274, decay_const_56Co, delta_decay_const_56Co, start_time_in_sec_Fe02, 2)
delta_A0_56Co_foil_3 = calc_delta_A0(1093.8536247846293, 1054.180513 , 75.49398, decay_const_56Co, delta_decay_const_56Co, start_time_in_sec_Fe03, 2)
delta_A0_56Co_foil_4 = calc_delta_A0(1394.464457609836, 1322.51407, 92.508566, decay_const_56Co, delta_decay_const_56Co, start_time_in_sec_Fe04, 2)
delta_A0_56Co_foil_5 = calc_delta_A0(2153.831008503944, 2041.734886, 141.361613, decay_const_56Co, delta_decay_const_56Co, start_time_in_sec_Fe05, 2)


print(f'56Co, Foil 1: delta_A0 = {delta_A0_56Co_foil_1}')
print(f'56Co, Foil 2: delta_A0 = {delta_A0_56Co_foil_2}')
print(f'56Co, Foil 3: delta_A0 = {delta_A0_56Co_foil_3}')
print(f'56Co, Foil 4: delta_A0 = {delta_A0_56Co_foil_4}')
print(f'56Co, Foil 5: delta_A0 = {delta_A0_56Co_foil_5}')

