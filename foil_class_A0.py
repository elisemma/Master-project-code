import numpy as np 
import matplotlib.pyplot as plt 
import curie as ci 
import pandas as pd 


class Foil_A0: 
    #wclass hich returns the beam current with unc for one foil
    def __init__(self, foil_name, reaction_product, A0, A0_unc):
        self.foil_name = foil_name
        self.reaction_product = reaction_product
        self.A0 = A0
        self.A0_unc = A0_unc
        # self.R = R 
        # self.R_unc = R_unc
        self.areal_dens = None
        self.areal_dens_unc_percent = None
        self.decay_const = None
        self.decay_const_unc = None
        # self.beam_energy_in_foil = None
        self.beam_current = None
        self.beam_current_unc = None
        self.weighted_average_beam_current = None
        self.var_weighted_average_beam_current = None
        self.calc_xs = None
        self.calc_xs_unc = None


    def assign_areal_dens_w_unc_percent(self): #[mg/cm2]
        areal_dens_dict = { 'Zr01': {'areal_dens': 16.142, 'areal_dens_unc_percent': 0.6381}, 
                            'Zr02': {'areal_dens': 16.378, 'areal_dens_unc_percent': 0.2992},
                            'Zr03': {'areal_dens': 16.140, 'areal_dens_unc_percent': 0.2478},
                            'Zr04': {'areal_dens': 16.195, 'areal_dens_unc_percent': 0.5619},
                            'Zr05': {'areal_dens': 16.400, 'areal_dens_unc_percent': 0.2195},
                            'Zr06': {'areal_dens': 16.177, 'areal_dens_unc_percent': 1.0571},
                            'Zr07': {'areal_dens': 16.009, 'areal_dens_unc_percent': 0.1124},
                            'Zr08': {'areal_dens': 16.246, 'areal_dens_unc_percent': 0.1662},
                            'Zr09': {'areal_dens': 16.225, 'areal_dens_unc_percent': 0.4253},
                            'Zr10': {'areal_dens': 16.237, 'areal_dens_unc_percent': 0.1478},

                            'Ti01': {'areal_dens': 11.535, 'areal_dens_unc_percent': 0.8236},
                            'Ti02': {'areal_dens': 11.740, 'areal_dens_unc_percent': 0.7751},
                            'Ti03': {'areal_dens': 11.273, 'areal_dens_unc_percent': 0.7629},
                            'Ti04': {'areal_dens': 11.099, 'areal_dens_unc_percent': 0.8289},
                            'Ti05': {'areal_dens': 11.317, 'areal_dens_unc_percent': 0.2651},
                            'Ti06': {'areal_dens': 10.990, 'areal_dens_unc_percent': 0.7734},
                            'Ti08': {'areal_dens': 11.169, 'areal_dens_unc_percent': 0.9319},
                            'Ti09': {'areal_dens': 11.234, 'areal_dens_unc_percent': 0.7655},
                            'Ti10': {'areal_dens': 11.239, 'areal_dens_unc_percent': 2.1443},
                            'Ti11': {'areal_dens': 11.001, 'areal_dens_unc_percent': 0.1364},

                            'Ni01': {'areal_dens': 23.253, 'areal_dens_unc_percent': 0.0645},
                            'Ni02': {'areal_dens': 23.103, 'areal_dens_unc_percent': 0.1688},
                            'Ni03': {'areal_dens': 23.064, 'areal_dens_unc_percent': 0.2992},
                            'Ni04': {'areal_dens': 23.201, 'areal_dens_unc_percent': 0.2974},
                            'Ni05': {'areal_dens': 22.746, 'areal_dens_unc_percent': 0.0835},

                            'Fe01': {'areal_dens': 20.061, 'areal_dens_unc_percent': 0.1495},
                            'Fe02': {'areal_dens': 20.061, 'areal_dens_unc_percent': 0.1538},
                            'Fe03': {'areal_dens': 20.187, 'areal_dens_unc_percent': 0.5548},
                            'Fe04': {'areal_dens': 20.100, 'areal_dens_unc_percent': 0.2488},
                            'Fe05': {'areal_dens': 20.108, 'areal_dens_unc_percent': 0.5570},}


        self.areal_dens = areal_dens_dict[self.foil_name]['areal_dens']
        self.areal_dens_unc_percent = areal_dens_dict[self.foil_name]['areal_dens_unc_percent']


    def assign_molar_mass(self):
        # Assigning the molar mass 
        molar_mass_dict = {'Zr': {'molar_mass': 91.2235, 'molar_mass_unc': 0.0005},
                           'Ni': {'molar_mass': 58.6934, 'molar_mass_unc': 0.0004},
                           'Ti': {'molar_mass': 47.867, 'molar_mass_unc': 0.001},
                           'Fe': {'molar_mass': 55.845, 'molar_mass_unc': 0.002}}

        self.molar_mass  = molar_mass_dict[self.foil_name[0:2]]['molar_mass']
        self.molar_mass_unc = molar_mass_dict[self.foil_name[0:2]]['molar_mass_unc']


    def calculate_decay_constant_w_unc(self):
        ip = ci.Isotope(self.reaction_product)
        half_life, half_life_unc = ip.half_life('s', True)
        # print(f't_1/2: {half_life} +- {half_life_unc}')

        self.decay_const = np.log(2)/half_life #1/[s])
        self.decay_const_unc = np.log(2)* half_life_unc /(half_life**2) #1/[s]


    def assign_beam_current_w_unc(self):
        # beam_current_dict = { 'Zr01': {'beam_current': 128.21543781323695, 'beam_current_unc': 12.74029303652092}, 
        #                       'Zr02': {'beam_current': 123.4575028703829, 'beam_current_unc': 7.710074388016526},
        #                       'Zr03': {'beam_current': 131.18779932476846, 'beam_current_unc': 4.0835044453116245},
        #                       'Zr04': {'beam_current': 96.63429803291486, 'beam_current_unc': 23.282870237984515},
        #                       'Zr05': {'beam_current': 7.443715398272095, 'beam_current_unc': 6.995006437642473},
        #                       'Ni01': {'beam_current': 136.2920634033693, 'beam_current_unc': 15.088128123068422},
        #                       'Ni02': {'beam_current': 129.94732270947284, 'beam_current_unc': 6.838926948100149},
        #                       'Ni03': {'beam_current': 134.93337063605304, 'beam_current_unc': 2.1254299196587896},
        #                       'Ni04': {'beam_current': 115.47046333281607, 'beam_current_unc': 5.3612072356381315},
        #                       'Ni05': {'beam_current': 7.443715398272095, 'beam_current_unc': 6.995006437642473},
        #                       'Ti01': {'beam_current': 121.3537108681938 , 'beam_current_unc': 2.074524514096691},
        #                       'Ti02': {'beam_current': 117.9039121079045, 'beam_current_unc': 1.8458020844093919},
        #                       'Ti03': {'beam_current': 127.3507935992085, 'beam_current_unc': 0.18696389527430857},
        #                       'Ti04': {'beam_current': 78.0621294846752, 'beam_current_unc': 18.801405804383073},
        #                       'Ti05': {'beam_current': 0, 'beam_current_unc': 0}}
        beam_current_dict = { 'Zr01': {'beam_current': 125.94248971622208, 'beam_current_unc': 3.9985510638858535},
                              'Zr02': {'beam_current': 129.10847496934917, 'beam_current_unc': 3.1486251610340217},
                              'Zr03': {'beam_current': 131.37180938909623, 'beam_current_unc': 2.097187106972314},
                              'Zr04': {'beam_current': 113.73832236122058, 'beam_current_unc': 5.0486614591625},
                              'Zr05': {'beam_current': 2.2556448798115443, 'beam_current_unc': 15.512092954286548},
                              'Zr06': {'beam_current': 90.40077287936163,  'beam_current_unc': 2.327839111192331},
                              'Zr07': {'beam_current': 90.93641924096387,  'beam_current_unc': 1.914904273219334},
                              'Zr08': {'beam_current': 94.82509814227382,  'beam_current_unc': 1.5491782774473961},
                              'Zr09': {'beam_current': 91.51235172180517,  'beam_current_unc': 2.44948286925908},
                              'Zr10': {'beam_current': 89.92515787653733,  'beam_current_unc': 2.1989658716372884},

                              'Ni01': {'beam_current': 133.93230394113667, 'beam_current_unc': 11.733851818232212},
                              'Ni02': {'beam_current': 126.03238410288209, 'beam_current_unc': 10.747827837439031},
                              'Ni03': {'beam_current': 131.63746302701534, 'beam_current_unc': 4.962623163414479},
                              'Ni04': {'beam_current': 110.71894968109665, 'beam_current_unc': 4.537838585108911},
                              'Ni05': {'beam_current': 0.6963785878966307, 'beam_current_unc': 2.834434564450914}, #Both 58Co and 61Cu, but this should not be used

                              'Ti01': {'beam_current': 121.40624351092747, 'beam_current_unc': 2.1995864938476406},
                              'Ti02': {'beam_current': 118.25453712836615, 'beam_current_unc': 1.716130752885589},
                              'Ti03': {'beam_current': 123.56651302304839, 'beam_current_unc': 2.3284985000035077},
                              'Ti04': {'beam_current': 765.51551559187187, 'beam_current_unc': 16.02184893105766},
                              'Ti05': {'beam_current': 108.57835015293583, 'beam_current_unc': 33.36232017668373},
                              'Ti06': {'beam_current': 96.60791103979946,  'beam_current_unc': 3.508512656191025},
                              'Ti08': {'beam_current': 93.07819175880893,  'beam_current_unc': 1.0660775372399667},
                              'Ti09': {'beam_current': 95.97313899866144,  'beam_current_unc': 1.0109998126930624},
                              'Ti10': {'beam_current': 92.89981297143379,  'beam_current_unc': 0.821473240005668},
                              'Ti11': {'beam_current': 90.22608281896498,  'beam_current_unc': 0.5246886893150025},

                              'Fe01': {'beam_current': 91.71868901291089,  'beam_current_unc': 17.891687001061534},
                              'Fe02': {'beam_current': 95.78064880345202,  'beam_current_unc': 17.586591409360715},
                              'Fe03': {'beam_current': 96.48332515195193,  'beam_current_unc': 16.633314233940332},
                              'Fe04': {'beam_current': 86.86262139072429,  'beam_current_unc': 13.01251633646137},
                              'Fe05': {'beam_current': 77.50784615270966,  'beam_current_unc': 9.55268221176066}}

        self.beam_current = beam_current_dict[self.foil_name]['beam_current']
        self.beam_current_unc = beam_current_dict[self.foil_name]['beam_current_unc']



    def calculate_xs_w_unc(self):

        N_A = 6.0221408e+23 
        t_irr = 1200 #[s]
        t_irr_unc = 3 #[s]
        
        areal_dens = float(self.areal_dens)/1000.0 # g/cm2
        molar_dens = areal_dens/self.molar_mass # mol/cm2

        N_T_per_cm2 = N_A*molar_dens #[nuclei/cm^2] when areal_dens is given in mg/cm^2

        N_T = N_T_per_cm2*1.0e4 #[nuclei/m^2]

        areal_dens_unc = areal_dens*self.areal_dens_unc_percent/100
        N_T_unc = N_T*np.sqrt((areal_dens_unc/areal_dens)**2 + (self.molar_mass_unc/self.molar_mass)**2) #[nuclei/m^2]

        beam_current_in_d_per_s = self.beam_current/(1.60217634e-19*1.0e9)
        beam_current_in_d_per_s_unc = self.beam_current_unc/(1.60217634e-19*1.0e9)

        xs = self.R/(beam_current_in_d_per_s * N_T)*1e+28*1e3 #[mb]

        self.calc_xs = xs

        xs_unc = xs * np.sqrt( (self.R_unc/self.R)**2 + (beam_current_in_d_per_s_unc/beam_current_in_d_per_s)**2 + (N_T_unc/N_T)**2 )#[mb]
        self.calc_xs_unc = xs_unc

        # print(self.foil_name, self.reaction_product)
        # print(f'R unc %: {self.R_unc/self.R*100}')
        # print(f'bc unc %: {beam_current_in_d_per_s_unc/beam_current_in_d_per_s*100}')
        # print(f'N_T unc %: {N_T_unc/N_T*100}')
        # print(f'xs: {xs} +- {xs_unc}\n')







    def calculate_xs_w_unc_old(self): #does not work with production rate as input in the class. Needs A0 w unc

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