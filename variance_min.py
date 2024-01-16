import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.interpolate import splev, splrep



class VarianceMin:
    def __init__(self, stack, beam_energy, projectile, target_material, reaction_list, A0, areal_dens, areal_dens_unc_percent):
        self.stack = stack
        self.beam_energy = beam_energy
        self.projectile = projectile
        self.reaction_list = reaction_list
        self.target_material = target_material
        self.A0 = A0
        self.areal_dens = areal_dens
        self.areal_dens_unc_percent = areal_dens_unc_percent
        self.molar_mass = None
        self.molar_mass_unc = None
        self.decay_const_list = [] 
        


    def assign_molar_mass(self):
        # Assigning the molar mass of the target material
        if self.target_material == 'Ni':
            self.molar_mass = 58.6934
            self.molar_mass_unc = 0.0004

        elif self.target_material == 'Ti':
            self.molar_mass = 47.867
            self.molar_mass_unc = 0.001

        else:
            print('Error: Unknown target material')


    def calculate_decay_const(self):
        # Calculating the decay constant of the production nuclei by checking the reaction happening 
        for reaction_product in self.reaction_list:
            if reaction_product == '46Sc':
                self.decay_const_list.append(np.log(2)/(83.79*24*3600)) #1/[s]

            elif reaction_product == '48V':
                self.decay_const_list.append(np.log(2)/(15.9735*24*3600)) #1/[s]

            elif reaction_product == '56Co':
                self.decay_const_list.append(np.log(2)/(77.236*24*3600)) #1/[s])

            elif reaction_product == '58Co':
                self.decay_const_list.append(np.log(2)/(70.86*24*3600)) #1/[s])

            elif reaction_product == '61Cu':
                self.decay_const_list.append(np.log(2)/(3.339*3600)) #1/[s])

            else:
                print('Error: Unknown monitor reaction')


    def energy_in_foil(self, N, E0, dp): 
        # Using the curie data package and the stack information to calculate the energy in all of the foils in the stack
        st = ci.Stack(self.stack, E0=E0, N=N, particle=self.projectile, dp = dp)
        st.saveas(f'stack_calc_E_{E0}_dp_{dp}.csv')
        st_pd = pd.read_csv(f'stack_calc_E_{E0}_dp_{dp}.csv')

        return st_pd


    def monitor_cross_section(self):
        # Using the IAEA data to get the monitor cross section for the energy in a foil
        filename = reaction_product + '.txt'
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

        E = np.array(E_list)
        xs = np.array(xs_list)
        unc_xs = np.array(xs_unc_list)
      
        weights = 1/unc_xs
      
        spl = splrep(E, xs, w=weights, k=4)

        return spl


    def calculate_beam_current_w_unc(self, Energy, decay_const):
        # Calculating the beam current with uncertainty by using the end of beam activity, areal density, molar mass, decay constant and monitor cross section
        N_A = 6.0221408e+23
        t_irr = 1200 #[s]
        t_irr_unc = 3 #[s]
        N_T = float(self.areal_dens)*N_A/self.molar_mass*10 #[nuclei/m^2] when areal_dens is given in mg/cm^2

        xs_mon = splev(Energy, self.monitor_cross_section())

        beam_current = self.A0/(N_T*xs_mon*(1-np.exp(-decay_const*t_irr))) #[d/s]
        beam_current_in_A = beam_current*1.60217634e-19 #[A]
        beam_current_in_nA = beam_current_in_A*1e9 #[nA]


        areal_dens_unc = self.areal_dens*10*self.areal_dens_unc_percent/100
        N_T_unc = N_T*np.sqrt((areal_dens_unc/self.areal_dens)**2 + (self.molar_mass_unc/self.molar_mass)**2) #[nuclei/cm^2]
    
        dfdx_list = [] #Jacobian
        unc_list = []
    
        dA0 = self.A0*1e-8
        dfdA0 = (self.A0+dA0/(N_T*xs_mon*(1-np.exp(-decay_const*t_irr))) - self.A0-dA0/(N_T*xs_mon*(1-np.exp( -decay_const*t_irr))))/dA0
        dfdx_list.append(dfdA0)
        unc_list.append(A0_unc)
    
        dN_T = N_T*1e-8
        dfdN_T = (self.A0/((N_T+dN_T)*xs_mon*(1-np.exp(-decay_const*t_irr))) - self.A0/((N_T-dN_T)*xs_mon*(1- np.exp(-decay_const*t_irr))))/dN_T
        dfdx_list.append(dfdN_T)
        unc_list.append(N_T_unc)
    
        dxs = xs_mon*1e-8
        dfdxs = (self.A0/(N_T*(xs_mon+dxs)*(1-np.exp(-decay_const*t_irr))) - self.A0/(N_T*(xs_mon-dxs)*(1-np. exp(-decay_const*t_irr))))/dxs
        dfdx_list.append(dfdxs)
        unc_list.append(xs_mon_unc)
    
        dt_irr = t_irr*1e-8
        dfdt_irr = (self.A0/(N_T*xs_mon*(1-np.exp(-decay_const*(t_irr+dt_irr)))) - self.A0/(N_T*xs_mon*(1-np. exp(-decay_const*(t_irr-dt_irr)))))/dt_irr 
        dfdx_list.append(dfdt_irr)
        unc_list.append(t_irr_unc)
    
        dfdx = np.array(dfdx_list)
        unc = np.array(unc_list)
        beam_current_unc = np.sqrt(np.sum(np.multiply(dfdx,dfdx)* np.multiply(unc,unc)))*1. 60217634e-19*1e9
    
        return beam_current_in_nA, beam_current_unc


    def calculate_chi_squared(self):
        # Finding the chi squared of the beam current by giving the calculated beam current with uncertainty and the weighted average beam current which will be assigned as the true value
        calculated_beam_currents, beam_currents_unc = self.calculate_beam_current_w_unc()
        weighted_average_beam_current = self.calculate_weighted_average_beam_current()

        weighted_average_beam_current_array = np.zeros(len(calculated_beam_currents))
        weighted_average_beam_current_array.fill(weighted_average_beam_current)
        x_diff = weighted_average_beam_current_array-calculated_beam_currents 

        chi_squared = np.sum( np.multiply(x_diff, x_diff)/ np.multiply(beam_currents_unc, beam_currents_unc))

        return chi_squared


    def calculate_weighted_average_beam_current(self):
        # Calculating the weighted average beam current of all the monitor reactions happening in a specific foil
        beam_current_list = []
        beam_current_unc_list = []
        weight_list = []
    
        for i in range(len(self.reaction_list)): 
            beam_current, beam_current_unc = self.calculate_beam_current_w_unc(Energy, self.decay_const_list[i])
            beam_current_list.append(beam_current)
            beam_current_unc_list.append(beam_current_unc)
            weight = 1/beam_current_unc**2
            weight_list.append(weight)
    
        weighted_average_beam_cur = np.sum(np.array(beam_current_list)*np.array(weight_list))/np.   sum(np.array(weight_list))
        var_weighted_average_beam_cur = np.sum((beam_current_list-weighted_average_beam_cur)**2)/len(beam_current_list)
    
        chi2 = self.calculate_chi_squared()

        return weighted_average_beam_cur, var_weighted_average_beam_cur, chi_squared, beam_current_list, beam_current_unc_list



    def find_optimal_beam_current(self):
        # Varying the beam energy and areal densities of the targets in the stack to get the optimal beam current by using the chi squared minimization method























