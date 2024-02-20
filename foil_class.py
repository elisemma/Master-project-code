import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 
# from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator

#One file for running the stack calulations with varying dp values and save the csv files in a folder

#This file with one class for foil and one function for variance minimization


class Foil: 
    #wclass hich returns the beam current with unc for one foil
    def __init__(self, foil_name, target_material, reaction_list, A0_list, A0_unc_list, areal_dens, dp, scaling_factor):
        self.foil_name = foil_name
        self.target_material = target_material
        self.reaction_list = reaction_list
        self.A0_list = A0_list
        self.A0_unc_list = A0_unc_list
        self.areal_dens = areal_dens
        self.areal_dens_unc_percent = None
        self.decay_const_list = []
        self.decay_const_unc_list = []
        self.xs_mon_list = []
        self.xs_mon_unc_list = []
        self.beam_energy_in_foil = None
        self.beam_current_list = []
        self.beam_current_unc_list = []
        self.weighted_average_beam_current = None
        self.var_weighted_average_beam_current = None
        self.dp = dp
        self.scaling_factor = scaling_factor
        self.calc_xs_list = []
        self.calc_xs_unc_list = []


    def assign_areal_dens_unc_percent(self):
        areal_dens_unc_dict = {'Ni01': 0.0645, 'Ni02': 0.1688, 'Ni03': 0.2992, 'Ni04': 0.2974, 'Ni05': 0.0835,
                               'Ti01': 0.8236, 'Ti02': 0.7751, 'Ti03': 0.7629, 'Ti04': 0.8289, 'Ti05': 0.2651}
        self.areal_dens_unc_percent = areal_dens_unc_dict[self.foil_name]

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
                self.decay_const_unc_list.append(np.log(2)*0.04*24*3600/((83.79*24*3600)**2)) #1/[s]

            elif reaction_product == '48V':
                self.decay_const_list.append(np.log(2)/(15.9735*24*3600)) #1/[s]
                self.decay_const_unc_list.append(np.log(2)*0.0025*24*3600/((15.9735*24*3600)**2)) #1/[s]

            elif reaction_product == '56CO':
                self.decay_const_list.append(np.log(2)/(77.236*24*3600)) #1/[s])
                self.decay_const_unc_list.append(np.log(2)*0.026*24*3600/((77.236*24*3600)**2)) #1/[s]

            elif reaction_product == '58CO':
                self.decay_const_list.append(np.log(2)/(70.86*24*3600)) #1/[s])
                self.decay_const_unc_list.append(np.log(2)*0.06*24*3600/((70.86*24*3600)**2)) #1/[s]

            elif reaction_product == '61CU':
                self.decay_const_list.append(np.log(2)/(3.339*3600)) #1/[s])
                self.decay_const_unc_list.append(np.log(2)*0.008*3600/((3.339*3600)**2)) #1/[s]

            else:
                print('Error: Unknown monitor reaction')



    def find_monitor_cross_section(self):
        #Importing the energy/fluxes from the stack calculation
        flux_file = f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/Stack_calculations/stack_30MeV_dp_{self.dp:.3f}_fluxes.csv'
        csv_flux_data = pd.read_csv(flux_file)
        
        target_flux_data = csv_flux_data.loc[csv_flux_data['name'] == self.foil_name]
        energy = target_flux_data.loc[:,'energy']
        flux = target_flux_data.loc[:,'flux']
        energy_from_stack_calc_list = energy.values.tolist()
        flux_from_stack_calc_list = flux.values.tolist() #[MeV]
        # energy_from_stack_calc = energy.values.tolist()
        # flux_from_stack_calc = flux.values.tolist() #[MeV]

        # Deleting the pileups near E=0
        energy_from_stack_calc_unfiltered = np.array(energy_from_stack_calc_list)
        flux_from_stack_calc_unfiltered = np.array(flux_from_stack_calc_list)
        energy_from_stack_calc = energy_from_stack_calc_unfiltered[energy_from_stack_calc_unfiltered >= 0.275]
        flux_from_stack_calc = flux_from_stack_calc_unfiltered[energy_from_stack_calc_unfiltered >= 0.275]

        self.beam_energy_in_foil = np.sum(energy_from_stack_calc * flux_from_stack_calc)/np.sum(flux_from_stack_calc)
        # Using the IAEA data to get the monitor cross section for the energy in a foil
        for reaction_product in self.reaction_list:
            filename = './Monitor_cross_section_data/IAEA_monitor_xs_' + reaction_product + '.txt'
            with open(filename) as file:
                lines = file.readlines()[6:-1]
                E_mon_list = []
                xs_list = []
                xs_unc_list = []

                for line in lines:
                    words = line.split()
                    E_mon_list.append(float(words[0]))
                    xs_list.append(float(words[1]))
                    xs_unc_list.append(float(words[2]))
                file.close()


            E_mon = np.array(E_mon_list)
            xs_mon = np.array(xs_list)*1e-28*1e-3 #convert mb to m^2
            unc_xs_mon = np.array(xs_unc_list)*1e-28*1e-3 #convert mb to m^2
            weights = 1 / unc_xs_mon
            # interp_xs = interp1d(E_mon, xs_mon,kind='linear')
            # interp_unc_xs = interp1d(E_mon, unc_xs_mon, kind='linear')
            interp_xs = PchipInterpolator(E_mon, xs_mon)
            interp_unc_xs = PchipInterpolator(E_mon, unc_xs_mon)


            # #Plotting to check theinterpolation by eye
            # energy_plotting_array= np.linspace(0,50,10000)
            # plt.plot(E_mon, xs_mon, 'bo', label='data')
            # plt.plot(energy_plotting_array, interp_xs(energy_plotting_array))
            # plt.legend()
            # plt.show()


            #Need to take the integral of the xs for the different energybins in the foil to find the right monitor cross section
            xs_for_energies_in_foil = interp_xs(energy_from_stack_calc)
            xs_un_norm = np.trapz(np.multiply(xs_for_energies_in_foil,flux_from_stack_calc), x=energy_from_stack_calc)
            xs_norm = xs_un_norm/np.trapz(flux_from_stack_calc, x= energy_from_stack_calc)
            self.xs_mon_list.append(xs_norm)

            #Finding the uncertainty in the same way as the cross section
            xs_unc_for_energies_in_foil = interp_unc_xs(energy_from_stack_calc)
            xs_unc_un_norm = np.trapz(np.multiply(xs_unc_for_energies_in_foil,flux_from_stack_calc), x=energy_from_stack_calc)
            xs_unc_norm = xs_unc_un_norm/np.trapz(flux_from_stack_calc, x= energy_from_stack_calc)
            self.xs_mon_unc_list.append(xs_unc_norm)





    def calculate_beam_currents_w_unc(self):
        # Calculating the beam current with uncertainty by using the end of beam activity, areal density, molar mass, decay constant and monitor cross section
        N_A = 6.0221408e+23
        t_irr = 1200 #[s]
        t_irr_unc = 3 #[s]
        # N_T_per_cm2 = float(self.areal_dens/1000)*N_A/self.molar_mass #[nuclei/cm^2] when areal_dens is given in mg/cm^2
        
        areal_dens = float(self.areal_dens)/1000.0 # g/cm2
        molar_dens = areal_dens/self.molar_mass # mol/cm2

        N_T_per_cm2 = N_A*molar_dens # nuclei / cm2 

        # N_T_per_cm2 = float(self.areal_dens/1000)*N_A/self.molar_mass #[nuclei/cm^2] when areal_dens is given in mg/cm^2
        N_T = N_T_per_cm2*1.0e4 #[nuclei/m^2]


        for i in range(len(self.reaction_list)):
            beam_current = self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*t_irr))) #[d/s]
            beam_current_in_A = beam_current*1.60217634e-19 #[A]
            beam_current_in_nA = beam_current_in_A*1.0e9 #[nA]

            areal_dens_unc = areal_dens*self.areal_dens_unc_percent/100
            N_T_unc = N_T*np.sqrt((areal_dens_unc/areal_dens)**2 + (self.molar_mass_unc/self.molar_mass)**2) #[nuclei/m^2]
        
            dfdx_list = [] #Jacobian
            unc_list = []
        
            dA0 = self.A0_list[i]*1e-8
            dfdA0 = (self.A0_list[i]+dA0/(N_T*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]-dA0/(N_T*self.xs_mon_list[i]*(1-np.exp( -self.decay_const_list[i]*t_irr))))/dA0
            dfdx_list.append(dfdA0)
            unc_list.append(self.A0_unc_list[i])
        
            dN_T = N_T*1e-8
            dfdN_T = (self.A0_list[i]/((N_T+dN_T)*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]/((N_T-dN_T)*self.xs_mon_list[i]*(1- np.exp(-self.decay_const_list[i]*t_irr))))/dN_T
            dfdx_list.append(dfdN_T)
            unc_list.append(N_T_unc)
        
            dxs = self.xs_mon_list[i]*1e-8
            dfdxs = (self.A0_list[i]/(N_T*(self.xs_mon_list[i]+dxs)*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]/(N_T*(self.xs_mon_list[i]-dxs)*(1-np. exp(-self.decay_const_list[i]*t_irr))))/dxs
            dfdx_list.append(dfdxs)
            unc_list.append(self.xs_mon_unc_list[i])
        
            dt_irr = t_irr*1e-8
            dfdt_irr = (self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np.exp(-self.decay_const_list[i]*(t_irr+dt_irr)))) - self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np. exp(-self.decay_const_list[i]*(t_irr-dt_irr)))))/dt_irr 
            dfdx_list.append(dfdt_irr)
            unc_list.append(t_irr_unc)

            ddecay_const = self.decay_const_list[i]*1e-8
            dfddecay_const = (self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np.exp(-(self.decay_const_list[i]+ddecay_const)*t_irr))) - self.A0_list[i]/(N_T*self.xs_mon_list[i]*(1-np. exp(-(self.decay_const_list[i]-ddecay_const)*t_irr))))/ddecay_const 
            dfdx_list.append(dfddecay_const)
            unc_list.append(self.decay_const_unc_list[i])
        
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



    def convert_beam_current_back_to_xs_w_unc(self):

        # Calculating the cross section with uncertainty by using the end of beam activity, areal density, molar mass, decay constant and beam current
        N_A = 6.0221408e+23
        t_irr = 1200 #[s]
        t_irr_unc = 3 #[s]
        # N_T_per_cm2 = float(self.areal_dens/1000)*N_A/self.molar_mass #[nuclei/cm^2] when areal_dens is given in mg/cm^2
        
        areal_dens = float(self.areal_dens)/1000.0 # g/cm2
        molar_dens = areal_dens/self.molar_mass # mol/cm2

        N_T_per_cm2 = N_A*molar_dens # nuclei / cm2 

        # N_T_per_cm2 = float(self.areal_dens/1000)*N_A/self.molar_mass #[nuclei/cm^2] when areal_dens is given in mg/cm^2
        N_T = N_T_per_cm2*1.0e4 #[nuclei/m^2]

        beam_current_in_d_per_s = self.weighted_average_beam_current/(1.60217634e-19*1.0e9)*self.scaling_factor
        beam_current_in_d_per_s_unc = np.sqrt(self.var_weighted_average_beam_current)/(1.60217634e-19*1.0e9)*self.scaling_factor



        for i in range(len(self.reaction_list)):
            
            xs = self.A0_list[i]/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const_list[i]*t_irr)))*1e+28*1e3 #[mb]

            areal_dens_unc = areal_dens*self.areal_dens_unc_percent/100
            N_T_unc = N_T*np.sqrt((areal_dens_unc/areal_dens)**2 + (self.molar_mass_unc/self.molar_mass)**2) #[nuclei/m^2]
        
            dfdx_list = [] #Jacobian
            unc_list = []
        
            dA0 = self.A0_list[i]*1e-8
            dfdA0 = (self.A0_list[i]+dA0/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]-dA0/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const_list[i]*t_irr))))/dA0
            dfdx_list.append(dfdA0)
            unc_list.append(self.A0_unc_list[i])
        
            dN_T = N_T*1e-8
            dfdN_T = (self.A0_list[i]/((N_T+dN_T)*beam_current_in_d_per_s*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]/((N_T-dN_T)*beam_current_in_d_per_s*(1- np.exp(-self.decay_const_list[i]*t_irr))))/dN_T
            dfdx_list.append(dfdN_T)
            unc_list.append(N_T_unc)
        
            dbc = self.beam_current_list[i]*1e-8
            dfdbc = (self.A0_list[i]/(N_T*(beam_current_in_d_per_s+dbc)*(1-np.exp(-self.decay_const_list[i]*t_irr))) - self.A0_list[i]/(N_T*(beam_current_in_d_per_s-dbc)*(1-np.exp(-self.decay_const_list[i]*t_irr))))/dbc
            dfdx_list.append(dfdbc)
            unc_list.append(self.beam_current_unc_list[i])
        
            dt_irr = t_irr*1e-8
            dfdt_irr = (self.A0_list[i]/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const_list[i]*(t_irr+dt_irr)))) - self.A0_list[i]/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const_list[i]*(t_irr-dt_irr)))))/dt_irr
            dfdx_list.append(dfdt_irr)
            unc_list.append(t_irr_unc)

            ddecay_const = self.decay_const_list[i]*1e-8
            dfddecay_const = (self.A0_list[i]/(N_T*beam_current_in_d_per_s*(1-np.exp(-(self.decay_const_list[i]+ddecay_const)*t_irr))) - self.A0_list[i]/(N_T*beam_current_in_d_per_s*(1-np.exp(-(self.decay_const_list[i]-ddecay_const)*t_irr))))/ddecay_const
            dfdx_list.append(dfddecay_const)
            unc_list.append(self.decay_const_unc_list[i])
        
            dfdx = np.array(dfdx_list)
            unc = np.array(unc_list)
            xs_unc = np.sqrt(np.sum(np.multiply(dfdx,dfdx)* np.multiply(unc,unc)))*1e+28*1e3 #[mb]

            self.calc_xs_list.append(xs)
            self.calc_xs_unc_list.append(xs_unc)

    



