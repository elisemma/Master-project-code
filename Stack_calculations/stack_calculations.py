import numpy as np 
import curie as ci 



#print(ci.COMPOUND_LIST): ['RbCl', 'Inconel', 'H2O', 'SS_316', 'NaCl', 'Kapton', 'Acrylic', 'Air', 'Silicone', 'SrCO3', 'Polyester', 'Polyethylene', 'HDPE']

stack_30MeV = [{'compound':'SS_316', 'name': 'SS1', 'ad':100.199},

               {'compound':'Kapton', 'name': 'Kapton1', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone1', 't':0.043}, 
               {'compound':'Ni', 'name': 'Ni01', 'ad':23.253}, 
               {'compound':'Silicone', 'name': 'Silicone2', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton2', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton3', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone3', 't':0.043},
               {'compound':'Zr', 'name': 'Zr01', 'ad':16.142},
               {'compound':'Silicone', 'name': 'Silicone4', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton4', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton5', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone5', 't':0.043},
               {'compound':'Ti', 'name': 'Ti01', 'ad':11.535},
               {'compound':'Silicone', 'name': 'Silicone6', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton6', 'ad':3.6128},

               {'compound':'Al', 'name': 'E1', 'ad':68.312},



               {'compound':'Al', 'name': 'E2', 'ad':68.245},

               {'compound':'Kapton', 'name': 'Kapton7', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone7', 't':0.043},
               {'compound':'Ni', 'name': 'Ni02', 'ad':23.103},
               {'compound':'Silicone', 'name': 'Silicone8', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton8', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton9', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone9', 't':0.043},
               {'compound':'Zr', 'name': 'Zr02', 'ad':16.378},
               {'compound':'Silicone', 'name': 'Silicone10', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton10', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton11', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone11', 't':0.043},
               {'compound':'Ti', 'name': 'Ti02', 'ad':11.740},
               {'compound':'Silicone', 'name': 'Silicone12', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton12', 'ad':3.6128},



               {'compound':'Al', 'name': 'E3', 'ad':68.349},

               {'compound':'Kapton', 'name': 'Kapton13', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone13', 't':0.043},
               {'compound':'Ni', 'name': 'Ni03', 'ad':23.064},
               {'compound':'Silicone', 'name': 'Silicone14', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton14', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton15', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone15', 't':0.043},
               {'compound':'Zr', 'name': 'Zr03', 'ad':16.140},
               {'compound':'Silicone', 'name': 'Silicone16', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton16', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton17', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone17', 't':0.043},
               {'compound':'Ti', 'name': 'Ti03', 'ad':11.273},
               {'compound':'Silicone', 'name': 'Silicone18', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton18', 'ad':3.6128},



               {'compound':'Kapton', 'name': 'Kapton19', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone19', 't':0.043},
               {'compound':'Ni', 'name': 'Ni04', 'ad':23.201},
               {'compound':'Silicone', 'name': 'Silicone20', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton20', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton21', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone21', 't':0.043},
               {'compound':'Zr', 'name': 'Zr04', 'ad':16.195},
               {'compound':'Silicone', 'name': 'Silicone22', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton22', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton23', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone23', 't':0.043},
               {'compound':'Ti', 'name': 'Ti04', 'ad':11.099},
               {'compound':'Silicone', 'name': 'Silicone24', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton24', 'ad':3.6128},



               {'compound':'Kapton', 'name': 'Kapton25', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone25', 't':0.043},
               {'compound':'Ni', 'name': 'Ni05', 'ad':22.746},
               {'compound':'Silicone', 'name': 'Silicone26', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton26', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton27', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone27', 't':0.043},
               {'compound':'Zr', 'name': 'Zr05', 'ad':16.400},
               {'compound':'Silicone', 'name': 'Silicone28', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton28', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton29', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone29', 't':0.043},
               {'compound':'Ti', 'name': 'Ti05', 'ad':11.317},
               {'compound':'Silicone', 'name': 'Silicone30', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton30', 'ad':3.6128},

               {'compound':'SS_316', 'name': 'SS2', 'ad':100.865 }]




dp_array = np.arange(0.8, 1.21, 0.01)
# dp_array2 = np.arange(0.97, 0.99, 0.001)
# dp_array1 = np.arange(0.99, 1.001, 0.001)
# dp_array2 = np.arange(0.94, 0.96, 0.001)
# dp_array = np.union1d(dp_array1, dp_array2)
print(dp_array)
dp_array_length = len(dp_array)
index = 0

for dp in dp_array:
    index += 1
    print('__________________________________________')
    print(f'Running stack calculation for dp = {dp:.3f}')

    st = ci.Stack(stack_30MeV, E0=30, dE0 = 0.45, N=1e6, particle='d', dp = dp)
    st.saveas(f'stack_30MeV_dp_{dp:.3f}.csv')

    percent_done = index/dp_array_length*100
    print(f'{percent_done:.3f}% of the calculation is done')
    print('__________________________________________\n')





