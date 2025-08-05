import matplotlib.pylab as plt


# # A list of all isotopes you want extracted for plotting
# products= ['38-Sr- 90', 
#            '39-Y - 86', '39-Y - 87', '39-Y - 88', '39-Y - 95', 
#            '40-Zr- 86', '40-Zr- 87', '40-Zr- 88', '40-Zr- 89', '40-Zr- 95',
#            '41-Nb- 86', '41-Nb- 87', '41-Nb- 88', '41-Nb- 89', '41-Nb- 90', '41-Nb- 95', '41-Nb- 96'
#            ]
# # Target isotope - needed if you have multiple sub-directories for each target isotope in an element
# element = 'Zr'
# target = '96Zr'


products= ['25-Mn- 52', '25-Mn- 54', '25-Mn- 60',
           '26-Fe- 52', '26-Fe- 60', 
           '27-Co- 52', '27-Co- 55', '27-Co- 56', '27-Co- 57', '27-Co- 58', '27-Co- 60',
           '28-Ni- 52', '28-Ni- 55', '28-Ni- 56', '28-Ni- 57', '28-Ni- 65',
           '29-Cu- 55', '29-Cu- 57', '29-Cu- 61', '29-Cu- 64'
           ]
element = 'Ni'
target = '64Ni'


# products= ['20-Ca- 47', '20-Ca- 48'
#            '21-Sc- 44', '21-Sc- 46', '21-Sc- 47', '21-Sc- 48', 
#            '23-V - 48'
#            ]
# element = 'Ti'
# target = '50Ti'


# products= ['23-V - 48',
#            '24-Cr- 48',
#            '25-Mn- 48', '25-Mn- 52', '25-Mn- 54', 
#            '26-Fe- 48', '26-Fe- 52',
#            '27-Co- 48', '27-Co- 52', '27-Co- 55', '27-Co- 56', '27-Co- 57', '27-Co- 58'
#            ]
# element = 'Fe'
# target = '58Fe'


with open('./'+element+'/'+target+'/defaults.lst', 'r') as f:

    # energies = [2.0]
    energies = []
    cross_sections = dict.fromkeys(products)
    reading = False

    i=1

    for line in f:
        # print(i)
        # if ' Incident energy  ' in line and '(CMS)' not in line :
        if ' incident energy  ' in line and '(CMS)' not in line and 'C.M.' not in line:
            # print(line)
            # current_energy = float(line[18:28].replace('D',"E"))
            current_energy = float(line[40:].replace('D',"E").replace("MeV (LAB)\n",""))
            # print(current_energy)
            if len(energies) == 0 or current_energy != energies[-1]:
                energies.append(current_energy)
        elif 'production cross section' in line and 'reac:' in line:
            isotope = line[2:11]
            # print(isotope)
            if isotope == '23-V - 47':
                print(line)
            if isotope in cross_sections.keys():
                if cross_sections[isotope] is None:
                    cross_sections[isotope] = []
                while len(cross_sections[isotope])<len(energies)-1:
                    cross_sections[isotope].append(0.0)
                cross_sections[isotope].append(float(line[38:49]))
        # i=i+1


# print(cross_sections)


for isotope in products:
    # print(isotope)
    if cross_sections[isotope] is None:
        print('nothing for ', isotope)
        continue
    #print('plotting ', isotope)
    try:
        # print(energies)
        plt.plot(energies,cross_sections[isotope])
        plt.xlabel('Proton Energy (MeV)')
        plt.ylabel('Cross Section (mb)')
        plt.title(isotope)
        plt.savefig('./'+element+'/'+target+'/xs_plots/'+isotope+"G_empire.png")
        # plt.show()
        plt.close()
    except ValueError:
        continue
    with open('./'+element+'/'+target+'/xs_plots/'+isotope.replace(" ","")+"G_empire.txt",'w') as f:
    # with open('./'+element+'/'+target+'/xs_plots/'+isotope.replace("-","")+"G_empire.txt",'w') as f:
        for i in range(len(energies)):
            f.write(str(energies[i])+"\t"+str(cross_sections[isotope][i])+"\n")
