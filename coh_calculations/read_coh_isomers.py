import matplotlib.pylab as plt


# A list of all isotopes you want extracted for plotting
# products= ['052-116Te', '052-117Te', '052-118Te', '052-119Te', '052-120Te', '052-121Te','052-122Te',
#            '051-115Sb', '051-116Sb', '051-117Sb', '051-118Sb', '051-119Sb', '051-120Sb', '051-121Sb', 
#            '050-112Sn', '050-113Sn', '050-114Sn', '050-115Sn', '050-116Sn', '050-117Sn', '050-118Sn', '050-119Sn', '050-120Sn' #, 
#            # etc etc
#            ]

products = ['039-086Y', '039-088Y']
# Target isotope - needed if you have multiple sub-directories for each target isotope in an element
target = '96Zr'

# Run 'cat out_coh_*_MeV.dat > out_coh_merged.dat' in the ./output directory to unify files for parsing
with open('./'+target+'/output/out_coh_merged.dat','r') as f:

    energies = []
    cross_sections = dict.fromkeys(products)
    reading = False

    for line in f:
        if 'sum' in line:
            reading = False
        elif 'STABLE AND LONG-LIVED STATE PRODUCTIONS' in line:
            reading = True
        elif '#..................................' in line:
            current_energy = float(line.split()[-1])
            if len(energies) == 0 or current_energy != energies[-1]:
                energies.append(current_energy)
        elif reading:
            isotope = line.split()[1]
            if isotope in cross_sections.keys():
                if cross_sections[isotope] is None:
                    cross_sections[isotope] = {'G':[], 'M':[]}
                ex = float(line.split()[2])
                if ex > 0:
                    key = 'M'
                else:
                    key = 'G'
                while len(cross_sections[isotope][key])<len(energies)-1:
                    cross_sections[isotope][key].append(0.0)

                cross_sections[isotope][key].append(float(line.split()[5]))

for isotope in products:
    for key in ['M','G']:
        if cross_sections[isotope][key] is None:
            print('nothing for ', isotope)
            continue
        #print('plotting ', isotope)
        try:
            # Pad 0.0 to all runs with missing channel data
            while len(cross_sections[isotope][key])<len(energies):
                cross_sections[isotope][key].append(0.0)
            plt.plot(energies,cross_sections[isotope][key])
            plt.xlabel('Proton Energy (MeV)')
            plt.ylabel('Cross Section (mb)')
            plt.title(isotope+' '+key)
            # plt.savefig(isotope+key+"_coh.png")
            # plt.savefig('./'+target+'/plots/'+isotope+key+"_coh.png")
            plt.show()
            plt.close()
        except ValueError:
            continue
        # with open('files/'+isotope+key+"_coh.txt",'w') as f:
        with open('./'+target+'/'+isotope+key+"_coh.txt",'w') as f:

            for i in range(len(energies)):
                f.write(str(energies[i])+"\t"+str(cross_sections[isotope][key][i])+"\n")
