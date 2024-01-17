import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.interpolate import splev, splrep
# from statistics import variance
# import os 
# import glob
# import seaborn as sns



stack_30MeV = [{'compound':'Fe', 'name': 'SS1', 'ad':100.199},#stainless steel compund ci.COMPOUND_LIST can be printed to get the stainless steel compound, add kapton and silicone  
               {'compound':'Ni', 'name': 'Ni01', 'ad':23.253}, #recalc ad to check that they are correct 
               {'compound':'Zr', 'name': 'Zr01', 'ad':16.142},
               {'compound':'Ti', 'name': 'Ti01', 'ad':11.535},
               {'compound':'Al', 'name': 'E1', 'ad':68.312},

               {'compound':'Al', 'name': 'E2', 'ad':68.245},
               {'compound':'Ni', 'name': 'Ni02', 'ad':23.103},
               {'compound':'Zr', 'name': 'Zr02', 'ad':16.378},
               {'compound':'Ti', 'name': 'Ti02', 'ad':11.739},

               {'compound':'Al', 'name': 'E3', 'ad':68.349},
               {'compound':'Ni', 'name': 'Ni03', 'ad':23.064},
               {'compound':'Zr', 'name': 'Zr03', 'ad':16.141},
               {'compound':'Ti', 'name': 'Ti03', 'ad':11.273},

               {'compound':'Ni', 'name': 'Ni04', 'ad':23.201},
               {'compound':'Zr', 'name': 'Zr04', 'ad':16.195},
               {'compound':'Ti', 'name': 'Ti04', 'ad':11.098},

               {'compound':'Ni', 'name': 'Ni05', 'ad':22.746},
               {'compound':'Zr', 'name': 'Zr05', 'ad':16.400},
               {'compound':'Ti', 'name': 'Ti05', 'ad':11.317},
               {'compound':'Fe', 'name': 'SS2', 'ad':100.865}]


st = ci.Stack(stack_30MeV, E0=30, dE0 = 0.3, N=1000000, particle='d', dp = 1)
# st.plot()
# print(st.stack)



def get_xs_from_exfor_files(filename):
    with open(filename) as file:
        lines = file.readlines()[7:-1]
        E_list = []
        xs_list = []
        xs_unc_list = []

        for line in lines:
            words = line.split()
            E_list.append(float(words[0]))
            xs_list.append(float(words[1]))
            xs_unc_list.append(float(words[2]))
        file.close()

    # return np.array(E_list), np.array(E_unc_list), np.array(xs_list), np.array(xs_unc_list)
    return E_list, xs_list, xs_unc_list


# def interpol_xs_unc(xs_unc_list, E_xs_list):

#     E_xs_list, index = np.unique(np.array(E_xs_list), return_index = True)
#     xs_unc_list = np.array(xs_unc_list)
#     xs_unc_list = xs_unc_list[index]

#     cs = CubicSpline(E_xs_list, xs_unc_list)
#     # cs = interp1d(E_xs_list, xs_unc_list, kind = 'nearest')
#     # plt.plot(E_xs_list, cs(E_xs_list))
#     # plt.xlabel('Energy (MeV)')
#     # plt.ylabel('xs (mb)')
#     # plt.show()
#     return cs


def interpol_xs_unc(Energies, cross_sections, unc_cross_sections):
    E = np.array(Energies)
    xs = np.array(cross_sections)
    unc_xs = np.array(unc_cross_sections)

    weights = 1/unc_xs

    spl = splrep(E, xs, w=weights, k=4)

    return spl



def convert_A0_to_beam_current_w_unc(A0, A0_unc, areal_dens, areal_dens_unc_percent, Mm, Mm_unc, decay_const, xs_mon, xs_mon_unc):
    N_A = 6.0221408e+23
    t_irr = 1200 #[s]
    t_irr_unc = 3 #[s]
    N_T = float(areal_dens)*N_A/Mm*10 #[nuclei/m^2] when areal_dens is given in mg/cm^2

    beam_current = A0/(N_T*xs_mon*(1-np.exp(-decay_const*t_irr))) #[d/s]
    beam_current_in_A = beam_current*1.60217634e-19 #[A]
    beam_current_in_nA = beam_current_in_A*1e9 #[nA]

    areal_dens_unc = areal_dens*10*areal_dens_unc_percent/100
    N_T_unc = N_T*np.sqrt((areal_dens_unc/areal_dens)**2 + (Mm_unc/Mm)**2) #[nuclei/cm^2]

    dfdx_list = [] #Jacobian
    unc_list = []

    dA0 = A0*1e-8
    dfdA0 = (A0+dA0/(N_T*xs_mon*(1-np.exp(-decay_const*t_irr))) - A0-dA0/(N_T*xs_mon*(1-np.exp(-decay_const*t_irr))))/dA0
    dfdx_list.append(dfdA0)
    unc_list.append(A0_unc)

    dN_T = N_T*1e-8
    dfdN_T = (A0/((N_T+dN_T)*xs_mon*(1-np.exp(-decay_const*t_irr))) - A0/((N_T-dN_T)*xs_mon*(1-np.exp(-decay_const*t_irr))))/dN_T
    dfdx_list.append(dfdN_T)
    unc_list.append(N_T_unc)

    dxs = xs_mon*1e-8
    dfdxs = (A0/(N_T*(xs_mon+dxs)*(1-np.exp(-decay_const*t_irr))) - A0/(N_T*(xs_mon-dxs)*(1-np.exp(-decay_const*t_irr))))/dxs
    dfdx_list.append(dfdxs)
    unc_list.append(xs_mon_unc)

    dt_irr = t_irr*1e-8
    dfdt_irr = (A0/(N_T*xs_mon*(1-np.exp(-decay_const*(t_irr+dt_irr)))) - A0/(N_T*xs_mon*(1-np.exp(-decay_const*(t_irr-dt_irr)))))/dt_irr 
    dfdx_list.append(dfdt_irr)
    unc_list.append(t_irr_unc)

    dfdx = np.array(dfdx_list)
    unc = np.array(unc_list)
    beam_current_unc = np.sqrt(np.sum(np.multiply(dfdx,dfdx)* np.multiply(unc,unc)))*1.60217634e-19*1e9

    return beam_current_in_nA, beam_current_unc



def weighted_average_beam_current(reaction_list, A0_list, A0_unc_list, areal_dens_list,areal_dens_unc_list,  Mm_list, Mm_unc_list, decay_const_list, xs_mon_list, xs_mon_unc_list):
    beam_current_list = []
    beam_current_unc_list = []
    weight_list = []

    for i in range(len(reaction_list)):
        beam_current, beam_current_unc = convert_A0_to_beam_current_w_unc(A0_list[i], areal_dens_list[i], Mm_list[i],    
        decay_const_list[i], xs_mon_list[i], A0_unc_list[i], areal_dens_unc_list[i], Mm_unc_list[i], xs_mon_unc_list[i])
        beam_current_list.append(beam_current)
        beam_current_unc_list.append(beam_current_unc)
        weight = 1/beam_current_unc**2
        weight_list.append(weight)
        print(f'Beam current for {reaction_list[i]}: {beam_current:.3f} +- {beam_current_unc:.3f}nA')
        print(f'cross section for {reaction_list[i]}: {xs_mon_list[i]}\n')


    weighted_average_beam_cur = np.sum(np.array(beam_current_list)*np.array(weight_list))/np.sum(np.array(weight_list))
    var_weighted_average_beam_cur = np.sum((beam_current_list-weighted_average_beam_cur)**2)/len(beam_current_list)

    chi_squared = chi_squared_func(beam_current_list, beam_current_unc_list, weighted_average_beam_cur)

    print(f'Weighted average beam current: {weighted_average_beam_cur:.3f} +- {np.sqrt(var_weighted_average_beam_cur):.3f} nA')

    return weighted_average_beam_cur, var_weighted_average_beam_cur, chi_squared, beam_current_list, beam_current_unc_list




def chi_squared_func(x_list, x_unc_list, x_true):

    if (type(x_true)==float):
        x_true_array = np.zeros(len(x_list))
        x_true_array.fill(xs_true)
    else:
        x_true_array = np.array(x_true)

    x_array = np.array(x_list)
    x_unc_array = np.array(x_unc_list)
    x_diff = x_true_array-x_array

    chi_squared = np.sum( np.multiply(x_diff, x_diff)/ np.multiply(x_unc_array, x_unc_array))
    print('Nå er vi inni chi_squared_func________________________________________________________________')
    print(f'x_unc_array = {x_unc_array}')
    print(f'x_array = {x_array}')
    print(f'x_true_array = {x_true_array}')
    print(f'x_diff = {x_diff}')
    print(f'chi squared = {chi_squared}')
    print('______________________________________________________________________________________________')

    return chi_squared



# tregner en func til å variere E og dp for å finne optimal beam cur.
# def finding_minimal_chi2():
    #må tenke på dette




#_____________________________________________Running the code____________________________________________________________

# file = './Monitor_cross_section_data/natNi_dx_61Cu.txt'
file = './Monitor_cross_section_data/natTi_dx_46Sc.txt'

E_list, xs_list, xs_unc_list = get_xs_from_exfor_files(file)
xs_interpol = interpol_xs_unc(E_list, xs_list, xs_unc_list)
E_array = np.linspace(0,80,1000000)

# plt.plot(E_array, xs_interpol(E_array))
# plt.plot(E_array, splev(E_array, xs_interpol))
# plt.plot(E_list, xs_list)
plt.show()


