import numpy as np 
from datetime import datetime
import curie as ci




def calc_decay_constant(t_half):
    return np.log(2)/t_half

def calc_delta_decay_const(t_half, delta_t_half):
    delta_decay_const = np.log(2)*delta_t_half/(t_half**2)
    return delta_decay_const

def calc_delta_A0(A0, A_t, delta_A_t, decay_const, delta_decay_const, t, delta_t):
    delta_A0 = A0*np.sqrt((delta_A_t/A_t)**2 + decay_const**2*t**2*((delta_decay_const/decay_const)**2+(delta_t/t)**2))
    return delta_A0



t_half_58Co = 70.86*24*3600 #[s]
delta_t_half_58Co = 0.06*24*3600 #[s]
delta_decay_const_58Co = calc_delta_decay_const(t_half_58Co, delta_t_half_58Co)
decay_const_58Co = calc_decay_constant(t_half_58Co)

t_half_61Cu = 3.339*3600 #[s]
delta_t_half_61Cu = 0.008*3600 #[s]
delta_decay_const_61Cu = calc_delta_decay_const(t_half_61Cu, delta_t_half_61Cu)
decay_const_61Cu = calc_decay_constant(t_half_61Cu)

t_half_46Sc = 83.79*24*3600 #[s]
delta_t_half_46Sc = 0.04*24*3600 #[s]
delta_decay_const_46Sc = calc_delta_decay_const(t_half_46Sc, delta_t_half_46Sc)
decay_const_46Sc = calc_decay_constant(t_half_46Sc)

ip_47Sc = ci.Isotope('47SC')
t_half_47Sc, delta_t_half_47Sc = ip_47Sc.half_life('s', True)
delta_decay_const_47Sc = calc_delta_decay_const(t_half_47Sc, delta_t_half_47Sc)
decay_const_47Sc = calc_decay_constant(t_half_47Sc)

time_EOB = datetime(2017, 2, 13, 14, 27, 00)
start_time_Ni01 = datetime(2017, 3, 20, 13, 23, 21) # 03/20/2017 13:23:21
start_time_Ni02 = datetime(2017, 3, 20, 14, 26, 41) # 03/20/2017 14:26:41
start_time_Ni04 = datetime(2017, 3, 23, 15, 58, 9) # 03/23/2017 15:58:09
start_time_BENi05 = datetime(2017, 2, 13, 21, 37, 39) # 02/13/2017 21:37:39
start_time_DGNi05 = datetime(2017, 3, 24, 11, 00, 16) # 03/24/2017 11:00:16

start_time_Ti02 = datetime(2017, 3, 6, 14, 41, 00) # 03/06/2017 14:41:00
start_time_Ti03 = datetime(2017, 3, 1, 16, 00, 4) # 03/01/2017 16:00:04
start_time_Ti05 = datetime(2017, 3, 2, 14, 50, 40) # 03/02/2017 14:50:40




start_time_in_sec_Ni01 = (start_time_Ni01 - time_EOB).total_seconds()
start_time_in_sec_Ni02 = (start_time_Ni02 - time_EOB).total_seconds()
start_time_in_sec_Ni04 = (start_time_Ni04 - time_EOB).total_seconds()
start_time_in_sec_BENi05 = (start_time_BENi05 - time_EOB).total_seconds()
start_time_in_sec_DGNi05 = (start_time_DGNi05 - time_EOB).total_seconds()
start_time_in_sec_Ti02 = (start_time_Ti02 - time_EOB).total_seconds()
start_time_in_sec_Ti03 = (start_time_Ti03 - time_EOB).total_seconds()
start_time_in_sec_Ti05 = (start_time_Ti05 - time_EOB).total_seconds()


#                      calc_delta_A0(A0, A_t, delta_A_t, decay_const, delta_decay_const, t, delta_t)
# delta_A0_58Co_foil_1 = calc_delta_A0(3736.2259003636364, 2653.641819, 146.630547, decay_const_58Co, delta_decay_const_58Co, start_time_in_sec_Ni01, 2)
# delta_A0_58Co_foil_2 = calc_delta_A0(5709.463064094319, 4053.244032, 590.375328, decay_const_58Co, delta_decay_const_58Co, start_time_in_sec_Ni02, 2)
# delta_A0_58Co_foil_4 = calc_delta_A0(1095.00048859536, 751.676675, 276.234471, decay_const_58Co, delta_decay_const_58Co, start_time_in_sec_Ni04, 2)
# delta_A0_58Co_foil_5 = calc_delta_A0(25.452633286057697,  17.384299, 1.683928, decay_const_58Co, delta_decay_const_58Co, start_time_in_sec_DGNi05, 2)
# delta_A0_61Cu_foil_5 = calc_delta_A0(461.7853246484494,  100.832212, 25.852977, decay_const_61Cu, delta_decay_const_61Cu, start_time_in_sec_BENi05, 2)
# delta_A0_46Sc_foil_5 = calc_delta_A0(1.7690322124431181,  1.531436, 0.329658, decay_const_46Sc, delta_decay_const_46Sc, start_time_in_sec_Ti05, 2)

delta_A0_47Sc_foil_2 = calc_delta_A0(1819.5642722877774,  23.293047, 3.294042, decay_const_47Sc, delta_decay_const_47Sc, start_time_in_sec_Ti02, 2)
delta_A0_47Sc_foil_3 = calc_delta_A0(628.4558452489525,  22.550619, 8.827815, decay_const_47Sc, delta_decay_const_47Sc, start_time_in_sec_Ti03, 2)
# delta_A0_47Sc_foil_5 = calc_delta_A0(29.875896900429762,  0.810585 , 0.193502, decay_const_47Sc, delta_decay_const_47Sc, start_time_in_sec_Ti05, 2)



# print(f'58Co, Foil 1: delta_A0 = {delta_A0_58Co_foil_1}')
# print(f'58Co, Foil 2: delta_A0 = {delta_A0_58Co_foil_2}')
# print(f'58Co, Foil 4: delta_A0 = {delta_A0_58Co_foil_4}')
# print(f'58Co, Foil 5: delta_A0 = {delta_A0_58Co_foil_5}')
# print(f'61Cu, Foil 5: delta_A0 = {delta_A0_61Cu_foil_5}')
# print(f'46Sc, Foil 5: delta_A0 = {delta_A0_46Sc_foil_5}')
print(f'47Sc, Foil 2: delta_A0 = {delta_A0_47Sc_foil_2}')
print(f'47Sc, Foil 3: delta_A0 = {delta_A0_47Sc_foil_3}')
# print(f'47Sc, Foil 5: delta_A0 = {delta_A0_47Sc_foil_5}')


