import numpy as np 
import curie as ci 

#print(ci.COMPOUND_LIST): ['RbCl', 'Inconel', 'H2O', 'SS_316', 'NaCl', 'Kapton', 'Acrylic', 'Air', 'Silicone', 'SrCO3', 'Polyester', 'Polyethylene', 'HDPE']

stack_50MeV = [{'compound':'SS_316', 'name': 'SS3', 'ad':100.476},

               {'compound':'Kapton', 'name': 'Kapton1', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone1', 't':0.043}, 
               {'compound':'Fe', 'name': 'Fe01', 'ad':20.061}, 
               {'compound':'Silicone', 'name': 'Silicone2', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton2', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton3', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone3', 't':0.043},
               {'compound':'Zr', 'name': 'Zr06', 'ad':16.177},
               {'compound':'Silicone', 'name': 'Silicone4', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton4', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton5', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone5', 't':0.043},
               {'compound':'Ti', 'name': 'Ti06', 'ad':10.990},
               {'compound':'Silicone', 'name': 'Silicone6', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton6', 'ad':3.6128},

               {'compound':'Al', 'name': 'C1', 'ad':261.480},


               {'compound':'Kapton', 'name': 'Kapton7', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone7', 't':0.043},
               {'compound':'Fe', 'name': 'Fe02', 'ad':20.161},
               {'compound':'Silicone', 'name': 'Silicone8', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton8', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton9', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone9', 't':0.043},
               {'compound':'Zr', 'name': 'Zr07', 'ad':16.009},
               {'compound':'Silicone', 'name': 'Silicone10', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton10', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton11', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone11', 't':0.043},
               {'compound':'Ti', 'name': 'Ti08', 'ad':11.160},
               {'compound':'Silicone', 'name': 'Silicone12', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton12', 'ad':3.6128},

               {'compound':'Al', 'name': 'E4', 'ad':68.290},
               {'compound':'Al', 'name': 'E5', 'ad':68.237},


               {'compound':'Kapton', 'name': 'Kapton13', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone13', 't':0.043},
               {'compound':'Fe', 'name': 'Fe03', 'ad':20.187},
               {'compound':'Silicone', 'name': 'Silicone14', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton14', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton15', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone15', 't':0.043},
               {'compound':'Zr', 'name': 'Zr08', 'ad':16.246},
               {'compound':'Silicone', 'name': 'Silicone16', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton16', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton17', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone17', 't':0.043},
               {'compound':'Ti', 'name': 'Ti09', 'ad':11.234},
               {'compound':'Silicone', 'name': 'Silicone18', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton18', 'ad':3.6128},

               {'compound':'Al', 'name': 'E6', 'ad':68.215},
               {'compound':'Al', 'name': 'E7', 'ad':68.185},


               {'compound':'Kapton', 'name': 'Kapton19', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone19', 't':0.043},
               {'compound':'Fe', 'name': 'Fe04', 'ad':20.100},
               {'compound':'Silicone', 'name': 'Silicone20', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton20', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton21', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone21', 't':0.043},
               {'compound':'Zr', 'name': 'Zr09', 'ad':16.225},
               {'compound':'Silicone', 'name': 'Silicone22', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton22', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton23', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone23', 't':0.043},
               {'compound':'Ti', 'name': 'Ti10', 'ad':11.239},
               {'compound':'Silicone', 'name': 'Silicone24', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton24', 'ad':3.6128},

               {'compound':'Al', 'name': 'E8', 'ad':68.155},
               {'compound':'Al', 'name': 'E9', 'ad':68.177},


               {'compound':'Kapton', 'name': 'Kapton25', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone25', 't':0.043},
               {'compound':'Fe', 'name': 'Fe05', 'ad':20.108},
               {'compound':'Silicone', 'name': 'Silicone26', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton26', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton27', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone27', 't':0.043},
               {'compound':'Zr', 'name': 'Zr10', 'ad':16.237},
               {'compound':'Silicone', 'name': 'Silicone28', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton28', 'ad':3.6128},

               {'compound':'Kapton', 'name': 'Kapton29', 'ad':3.6128},
               {'compound':'Silicone', 'name': 'Silicone29', 't':0.043},
               {'compound':'Ti', 'name': 'Ti11', 'ad':11.001},
               {'compound':'Silicone', 'name': 'Silicone30', 't':0.043},
               {'compound':'Kapton', 'name': 'Kapton30', 'ad':3.6128},

               {'compound':'SS_316', 'name': 'SS6', 'ad':100.992 }]




dp_array = np.arange(0.8, 1.21, 0.01)
print(dp_array)
dp_array_length = len(dp_array)
index = 0

for dp in dp_array:
    index += 1
    print('__________________________________________')
    print(f'Running stack calculation for dp = {dp:.2f}')
    st = ci.Stack(stack_50MeV, E0=50, dE0 = 0.75, N=1e6, particle='d', dp = dp) #15% dE0
    st.saveas(f'stack_50MeV_dp_{dp:.2f}.csv')
    percent_done = index/dp_array_length*100
    
    print(f'{percent_done:.2f}% of the calculation is done')
    print('__________________________________________\n')

# print(st.stack)
# st.plot()





