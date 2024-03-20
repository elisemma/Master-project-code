from __future__ import print_function
from __future__ import with_statement

import matplotlib.pylab as plt
import os
import glob


def plot_levels(filenames):
    """Function to plot cumulative levels given in a zvd file"""



    

    for files in filenames:
        experimental_energies = []
        experimental_levels = []
        theoretical_energies = {}
        theoretical_levels = {}
        reading_experiment = False
        
        with open(files, 'r') as f:
            for line in f:
                #print(reading_data, line)
                if line[0] == '#':
                    continue
                elif 'fun:' in line:
                    if 'Exp' in line and len(experimental_energies) == 0:
                        reading_experiment = True
                        nuc =line.split()[-1]
                    elif 'Integral' in line:
                        reading_experiment = False
                        theo_label = line.split()[-1]
                        if theo_label not in theoretical_energies.keys():
                            theoretical_levels[theo_label] = []
                            theoretical_energies[theo_label] = []
#                elif '//' in line:
#                    print('changing reading data from ',reading_data)
#                    reading_data *= -1
#                    print('to ',reading_data)

                else:
                    try:
                        line = line.split()
                        energy = float(line[0])/10**6
                        levels = float(line[1])
                    except ValueError:
                        continue
                    if reading_experiment:
                        experimental_energies.append(energy)
                        experimental_levels.append(levels)
                    else:
                        theoretical_energies[theo_label].append(energy)
                        theoretical_levels[theo_label].append(levels)

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
        plt.xlabel('Excitiation Energy [MeV]')
        plt.ylabel('Cumulative Levels')
        ax.set_yscale('log')
        plt.title('Cumul_Levels_'+nuc)
        plt.plot(experimental_energies,experimental_levels,label='Experimental')
        for key in theoretical_energies.keys():
            plt.plot(theoretical_energies[key],theoretical_levels[key],label=key)
        plt.legend()
        plt.savefig("./levden_plots/"+nuc+"_levels.png")
        #plt.show()

	
#if name == '__main__':
label='defaults'
model='EGSM'
filenames = glob.glob(label+'-NL_'+model+'_*.zvd')
print(filenames)
# filenames = [label+'-NL_'+model+'_Te119.zvd',label+'-NL_'+model+'_Te120.zvd',label+'-NL_'+model+'_Te121.zvd',
#               label+'-NL_'+model+'_Te122.zvd',label+'-NL_'+model+'_Sn115.zvd',label+'-NL_'+model+'_Sn116.zvd',
#               label+'-NL_'+model+'_Sn117.zvd',label+'-NL_'+model+'_Sn118.zvd',label+'-NL_'+model+'_Sn119.zvd',
# 			label+'-NL_'+model+'_Sn120.zvd']

plot_levels(filenames)










