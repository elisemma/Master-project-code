import matplotlib.pylab as plt


# A list of all isotopes you want extracted for plotting
products= ['52-Te-116', '52-Te-117', '52-Te-118', '52-Te-119', '52-Te-120', '52-Te-121', '52-Te-122', 
           '51-Sb-115', '51-Sb-116', '51-Sb-117', '51-Sb-118', '51-Sb-119', '51-Sb-120', '51-Sb-121', 
           '50-Sn-112', '50-Sn-113', '50-Sn-114', '50-Sn-115', '50-Sn-116', '50-Sn-117', '50-Sn-118', '50-Sn-119', '50-Sn-120' #, 
           # etc etc
           ]

# Target isotope - needed if you have multiple sub-directories for each target isotope in an element
target = '123Sb'

with open('./'+target+'/defaults.lst', 'r') as f:

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
        plt.savefig(target+'/plots/'+isotope+"_empire.png")
        plt.close()
    except ValueError:
        continue
    # with open('.plots/'+isotope.replace(" ","")+"_empire.txt",'w') as f:
    with open('./'+target+'/plots/'+isotope.replace("-","")+"_empire.txt",'w') as f:
        for i in range(len(energies)):
            f.write(str(energies[i])+"\t"+str(cross_sections[isotope][i])+"\n")
