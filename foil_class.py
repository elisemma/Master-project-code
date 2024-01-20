import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
from scipy.interpolate import splev, splrep

#One file for running the stack calulations with varying dp values and save the csv files in a folder

#This file with one class for foil and one function for variance minimization


class Foil: 
    #wclass hich returns the beam current with unc for one foil
    def __init__(self, beam_energy_in_foil, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, areal_dens_unc_percent):
        self.beam_energy_in_foil = beam_energy_in_foil
        self.target_material = target_material
        self.reaction_list = reaction_list
        self.A0_list = A0_list
        self.A0_unc_list = A0_unc_list
        self.areal_dens = areal_dens
        self.areal_dens_unc_percent = areal_dens_unc_percent
        self.decay_const_list = []
        self.xs_mon_list = []
        self.xs_mon_unc_list = [0,0,0]
        self.beam_current_list = []
        self.beam_current_unc_list = []
        self.weighted_average_beam_current = None
        self.var_weighted_average_beam_current = None


    def assign_molar_mass(self):
        # Assigning the molar mass of the target material
        if self.target_material == 'Ni':
            self.molar_mass  = 58.6934
            self.molar_mass_unc = 0.0004

        elif self.target_material == 'Ti':
            self.molar_mass = 47.867
            self.molar_mass_unc = 0.001

        else:
            print('Error: Unknown target material')


    def calculate_decay_constant(self):
        # Calculating the decay constant of the production nuclei by checking the reaction happening 
        for reaction_product in self.reaction_list:
            if reaction_product == '46SC':
                self.decay_const_list.append(np.log(2)/(83.79*24*3600)) #1/[s]

            elif reaction_product == '48V':
                self.decay_const_list.append(np.log(2)/(15.9735*24*3600)) #1/[s]

            elif reaction_product == '56CO':
                self.decay_const_list.append(np.log(2)/(77.236*24*3600)) #1/[s])

            elif reaction_product == '58CO':
                self.decay_const_list.append(np.log(2)/(70.86*24*3600)) #1/[s])

            elif reaction_product == '61CU':
                self.decay_const_list.append(np.log(2)/(3.339*3600)) #1/[s])

            else:
                print('Error: Unknown monitor reaction')



    def find_monitor_cross_section(self):
        # Using the IAEA data to get the monitor cross section for the energy in a foil
        for reaction_product in self.reaction_list:
            filename = './Monitor_cross_section_data/IAEA_monitor_xs_' + reaction_product + '.txt'
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
            weights = 1 / unc_xs
            spl = splrep(E, xs, w=weights, k=4)

            # Interpolate the entire list of energies
            interpolated_xs = splev(self.beam_energy_in_foil, spl)

            # Append the interpolated value to the xs_mon_list
            self.xs_mon_list.append(interpolated_xs)
            self.xs_mon_unc_list.append(np.interp(self.beam_energy_in_foil, E, unc_xs))  # Assuming linear interpolation for uncertainty

    def calculate_beam_currents_w_unc(self):
        # Calculating the beam current with uncertainty by using the end of beam activity, areal density, molar mass, decay constant and monitor cross section
        N_A = 6.0221408e+23
        t_irr = 1200 #[s]
        t_irr_unc = 3 #[s]
        N_T_per_cm2 = float(self.areal_dens)/1000*N_A/self.molar_mass*10 #[nuclei/cm^2] when areal_dens is given in g/cm^2
        N_T = N_T_per_cm2*1e4 #[nuclei/m^2]

        for i in range(len(self.reaction_list)):
            beam_current = self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*t_irr))) #[d/s]
            beam_current_in_A = beam_current*1.60217634e-19 #[A]
            beam_current_in_nA = beam_current_in_A*1e9 #[nA]


            areal_dens_unc = self.areal_dens*10*self.areal_dens_unc_percent/100
            N_T_unc = N_T*np.sqrt((areal_dens_unc/self.areal_dens)**2 + (self.molar_mass_unc/self.molar_mass)**2) #[nuclei/cm^2]
        
            dfdx_list = [] #Jacobian
            unc_list = []
        
            dA0 = self.A0_list[i]*1e-8
            dfdA0 = (self.A0_list[i]+dA0/(N_T*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]-dA0/(N_T*self.xs_mon_list[i]*(1-np.exp( -self.decay_const_list[i]*t_irr))))/dA0
            dfdx_list.append(dfdA0)
            unc_list.append(self.A0_unc_list[i])
        
            # dN_T = N_T*1e-8
            # dfdN_T = (self.A0_list[i]/((N_T+dN_T)*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]/((N_T-dN_T)*self.xs_mon_list[i]*(1- np.exp(-self.decay_const_list[i]*t_irr))))/dN_T
            # dfdx_list.append(dfdN_T)
            # unc_list.append(N_T_unc)
        
            dxs = self.xs_mon_list[i]*1e-8
            dfdxs = (self.A0_list[i]/(N_T*(self.xs_mon_list[i]+dxs)*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]/(N_T*(self.xs_mon_list[i]-dxs)*(1-np. exp(-self.decay_const_list[i]*t_irr))))/dxs
            dfdx_list.append(dfdxs)
            unc_list.append(self.xs_mon_unc_list[i])
        
            dt_irr = t_irr*1e-8
            dfdt_irr = (self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*(t_irr+dt_irr)))) - self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np. exp(-self.decay_const_list[i]*(t_irr-dt_irr)))))/dt_irr 
            dfdx_list.append(dfdt_irr)
            unc_list.append(t_irr_unc)
        
            dfdx = np.array(dfdx_list)
            unc = np.array(unc_list)
            beam_current_unc = np.sqrt(np.sum(np.multiply(dfdx,dfdx)* np.multiply(unc,unc)))*1.60217634e-19*1e9

            self.beam_current_list.append(beam_current_in_nA)
            self.beam_current_unc_list.append(beam_current_unc)
    


    def calculate_weighted_average_beam_current(self):
        # Calculating the weighted average beam current of all the monitor reactions happening in a specific foil
        weight_list = []

        for i in range(len(self.reaction_list)):
            weight = 1 / self.beam_current_unc_list[i]**2
            weight_list.append(weight)

        # Calculate the weighted average beam current
        self.weighted_average_beam_current = np.average(self.beam_current_list, weights=weight_list)

        # Calculate the variance using the weighted_average_beam_current
        variance = np.average((np.array(self.beam_current_list) - self.weighted_average_beam_current)**2, weights=weight_list)
        self.var_weighted_average_beam_current = float(variance)



