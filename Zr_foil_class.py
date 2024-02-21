import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 


class Zr_foil: 
    #wclass hich returns the beam current with unc for one foil
    def __init__(self, foil_name, reaction_product, A0, A0_unc):
        self.foil_name = foil_name
        self.reaction_product = reaction_product
        self.A0 = A0
        self.A0_unc = A0_unc
        self.areal_dens = None
        self.areal_dens_unc_percent = None
        self.decay_const = None
        self.decay_const_unc = None
        self.beam_energy_in_foil = None
        self.beam_current = None
        self.beam_current_unc = None
        self.weighted_average_beam_current = None
        self.var_weighted_average_beam_current = None
        self.calc_xs = None
        self.calc_xs_unc = None


    def assign_areal_dens_w_unc_percent(self):
        areal_dens_dict = { 'Zr01': {'areal_dens': 16.142, 'areal_dens_unc_percent': 0.6381},
                            'Zr02': {'areal_dens': 16.378, 'areal_dens_unc_percent': 0.2992},
                            'Zr03': {'areal_dens': 16.140, 'areal_dens_unc_percent': 0.2478},
                            'Zr04': {'areal_dens': 16.195, 'areal_dens_unc_percent': 0.5619},
                            'Zr05': {'areal_dens': 16.400, 'areal_dens_unc_percent': 0.2195}}

        self.areal_dens = areal_dens_dict[self.foil_name]['areal_dens']
        self.areal_dens_unc_percent = areal_dens_dict[self.foil_name]['areal_dens_unc_percent']


    def assign_molar_mass(self):
        # Assigning the molar mass 
        self.molar_mass  = 91.2235
        self.molar_mass_unc = 0.0005


    def calculate_decay_constant_w_unc(self):
        # Calculating the decay constant of the production nuclei by checking the reaction happening 
        half_life_dict = {'96NB': {'half_life': 23.35*3600, 'half_life_unc': 5*3600},
                          '87Y':  {'half_life': 79.8*3600, 'half_life_unc': 3*3600}}

        self.decay_const = np.log(2)/half_life_dict[self.reaction_product]['half_life'] #1/[s])
        self.decay_const_unc = np.log(2)* half_life_dict[self.reaction_product]['half_life_unc'] /(half_life_dict[self.reaction_product]['half_life']**2) #1/[s]


    def assign_beam_current_w_unc(self):
        beam_current_dict = { 'Zr01': {'beam_current': 128.21543781323695, 'beam_current_unc': 12.74029303652092},
                              'Zr02': {'beam_current': 123.4575028703829, 'beam_current_unc': 7.710074388016526},
                              'Zr03': {'beam_current': 131.18779932476846, 'beam_current_unc': 4.0835044453116245},
                              'Zr04': {'beam_current': 96.63429803291486, 'beam_current_unc': 23.282870237984515},
                              'Zr05': {'beam_current': 0, 'beam_current_unc': 0}}

        self.beam_current = beam_current_dict[self.foil_name]['beam_current']
        self.beam_current_unc = beam_current_dict[self.foil_name]['beam_current_unc']


    def calculate_xs_w_unc(self):

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

        beam_current_in_d_per_s = self.beam_current/(1.60217634e-19*1.0e9)
        beam_current_in_d_per_s_unc = self.beam_current_unc/(1.60217634e-19*1.0e9)

            

        xs = self.A0/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const*t_irr)))*1e+28*1e3 #[mb]
        self.calc_xs = xs



        areal_dens_unc = areal_dens*self.areal_dens_unc_percent/100
        N_T_unc = N_T*np.sqrt((areal_dens_unc/areal_dens)**2 + (self.molar_mass_unc/self.molar_mass)**2) #[nuclei/m^2]
    
        dfdx_list = [] #Jacobian
        unc_list = []
    
        dA0 = self.A0*1e-8
        dfdA0 = (self.A0+dA0/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const*t_irr))) - self.A0-dA0/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const*t_irr))))/dA0
        dfdx_list.append(dfdA0)
        unc_list.append(self.A0_unc)
    
        dN_T = N_T*1e-8
        dfdN_T = (self.A0/((N_T+dN_T)*beam_current_in_d_per_s*(1-np.exp(-self.decay_const*t_irr))) - self.A0/((N_T-dN_T)*beam_current_in_d_per_s*(1- np.exp(-self.decay_const*t_irr))))/dN_T
        dfdx_list.append(dfdN_T)
        unc_list.append(N_T_unc)
    
        dbc = beam_current_in_d_per_s*1e-8
        dfdbc = (self.A0/(N_T*(beam_current_in_d_per_s+dbc)*(1-np.exp(-self.decay_const*t_irr))) - self.A0/(N_T*(beam_current_in_d_per_s-dbc)*(1-np.exp(-self.decay_const*t_irr))))/dbc
        dfdx_list.append(dfdbc)
        unc_list.append(beam_current_in_d_per_s_unc)
    
        dt_irr = t_irr*1e-8
        dfdt_irr = (self.A0/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const*(t_irr+dt_irr)))) - self.A0/(N_T*beam_current_in_d_per_s*(1-np.exp(-self.decay_const*(t_irr-dt_irr)))))/dt_irr
        dfdx_list.append(dfdt_irr)
        unc_list.append(t_irr_unc)

        ddecay_const = self.decay_const*1e-8
        dfddecay_const = (self.A0/(N_T*beam_current_in_d_per_s*(1-np.exp(-(self.decay_const+ddecay_const)*t_irr))) - self.A0/(N_T*beam_current_in_d_per_s*(1-np.exp(-(self.decay_const-ddecay_const)*t_irr))))/ddecay_const
        dfdx_list.append(dfddecay_const)
        unc_list.append(self.decay_const_unc)
    
        dfdx = np.array(dfdx_list)
        unc = np.array(unc_list)
        xs_unc = np.sqrt(np.sum(np.multiply(dfdx,dfdx)* np.multiply(unc,unc)))*1e+28*1e3 #[mb]

        self.calc_xs_unc = xs_unc