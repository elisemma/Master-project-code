import matplotlib.pylab as plt


# For this, you'll need to look up the excitation energy for each of your isomers in Empire, 
# found in your defaults.lst file.   Try doing Control+F to look for "is an isomer", and grab the Ex
# once you find it for a given isomer
# 
# Syntax:  'Ex': 'isotope'
# products= {'0.0895':'50-Sn-119', '0.1636':'51-Sb-122'}

products= {'0.2183':'39-Y-86', '0.3808':'39-Y-87', ' 0.5878':'40-Zr-89', '0.6817':'39-Y-90', '0.1355':'41-Nb-92', '0.2357':'41-Nb-95'}
element = 'Zr'
target = '96Zr'

# products= {'0.2712':'21-Sc-44'}
# element = 'Ti'
# target = '50Ti'

# products= {'0.3777':'25-Mn-52'}
# element = 'Fe'
# target = '57Fe'

# products= {'0.3777':'25-Mn-52', '0.0586':'27-Co-60'}
# element = 'Ni'
# target = '64Ni'

# # Target isotope - needed if you have multiple sub-directories for each target isotope in an element
# target = '123Sb'
# with open('./'+target+'/defaults.lst', 'r') as f:
# Target isotope - needed if you have multiple sub-directories for each target isotope in an element


with open('./'+element+'/'+target+'/defaults.lst', 'r') as f:

    # energies = [2.0]
    energies = []
    cross_sections = dict.fromkeys(products.keys())
    reading = False

    for line in f:
        if ' incident energy  ' in line and '(CMS)' not in line and 'C.M.' not in line :
            # print(line)
            # current_energy = float(line[18:28].replace('D',"E"))
            current_energy = float(line[40:].replace('D',"E").replace("MeV (LAB)\n",""))
            # print(current_energy)
            if len(energies) == 0 or current_energy != energies[-1]:
                energies.append(current_energy)
        elif 'is an isomer' in line and 'Level of energy' in line:
            energy = line.split()[3]
            #print(energy)
            if energy in cross_sections.keys():
                # print(line)
                if cross_sections[energy] is None:
                    cross_sections[energy] = []
                while len(cross_sections[energy])<len(energies)-1:
                    cross_sections[energy].append(0.0)
                cross_sections[energy].append(float(line.split()[11]))





for isotope in products:
    if cross_sections[isotope] is None:
        print('nothing for ', isotope)
        continue
    #print('plotting ', isotope)
    # print(energies)
    # print(cross_sections[isotope])
    try:
        plt.plot(energies,cross_sections[isotope])
        plt.xlabel('Proton Energy (MeV)')
        plt.ylabel('Cross Section (mb)')
        plt.title(products[isotope]+" M")
        plt.savefig('./'+element+'/'+target+'/xs_plots/'+products[isotope]+"M_empire.png")
        # plt.show()
        plt.close()
    except ValueError:
        continue
    # with open('./'+target+'/plots/'+products[isotope].replace(" ","")+"M_empire.txt",'w') as f:
    with open('./'+element+'/'+target+'/xs_plots/'+products[isotope].replace(" ","")+"M_empire.txt",'w') as f:
        for i in range(len(energies)):
            f.write(str(energies[i])+"\t"+str(cross_sections[isotope][i])+"\n")
