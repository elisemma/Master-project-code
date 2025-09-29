import numpy as np
import matplotlib.pyplot as plt




#86Sr(d,2n)86gY(cum)
E_86Sr_d_2n_86Y_Uddin = np.array([0, 10, 10.08, 10.47414, 10.86828, 11.26242, 11.65657, 12.05071, 12.44485, 12.83899, 13.23313, 13.62727, 14.02141, 14.41556, 14.8097, 15.20384, 15.59798, 15.99212, 16.38626, 16.7804, 17.17455, 17.56869, 17.96283, 18.35697, 18.75111, 19.14525, 19.53939, 19.93354, 20.32768, 20.72182, 21.11596, 21.5101, 21.90424, 22.29838, 22.69253, 23.08667, 23.48081, 23.87495, 24.26909, 24.66323, 25.05737, 25.45152, 25.84566, 26.2398, 26.63394, 27.02808, 27.42222, 27.81636, 28.21051, 28.60465, 28.99879, 29.39293, 29.78707, 30.18121, 30.57535, 30.96949, 31.36364, 31.75778, 32.15192, 32.54606, 32.9402, 33.33434, 33.72848]) #[MeV]
Yield_86Sr_d_2n_86Y_Uddin = np.array([0, 0.002099527, 0.484959413, 3.543875229, 8.727507359, 15.90074458, 24.92586673, 35.66380824, 47.97556036, 61.72275359, 76.76862422, 92.97876979, 110.2220372, 128.370474, 147.3002101, 166.8924321, 187.0330545, 207.6132698, 228.5298147, 249.685452, 270.9883407, 292.3526639, 313.6992792, 334.9549409, 356.0525744, 376.9312498, 397.5363844, 417.818874, 437.7355002, 457.2493113, 476.3286892, 494.9474267, 513.0845311, 530.7242341, 547.8550896, 564.470191, 580.5673423, 596.1481424, 611.2179321, 625.7855137, 639.8630464, 653.4652165, 666.6093179, 679.3152879, 691.6048987, 703.5016351, 715.0304017, 726.2173769, 737.0892989, 747.67348, 757.9977985, 768.0900209, 777.9776822, 787.6878374, 797.2468229, 806.6801495, 816.0119307, 825.2649194, 834.4605594, 843.6184571, 852.7563457, 861.8899452, 871.0329591]) #[MBq/uA h]


#86Sr(p,n)86gY(cum)
E_86Sr_p_n_86Y_Uddin = np.array([0, 5.66606, 6.12686, 6.58765, 7.04844, 7.50924, 7.97003, 8.43083, 8.89162, 9.35242, 9.81321, 10.27401, 10.7348, 11.1956, 11.65639, 12.11719, 12.57798, 13.03877, 13.49957, 13.96036, 14.42116, 14.88195, 15.34275, 15.80354, 16.26434, 16.72513, 17.18593, 17.64672, 18.10751, 18.56831, 19.0291, 19.4899, 19.95069, 20.41149, 20.87228, 21.33308, 21.79387, 22.25467, 22.71546, 23.17626, 23.63705, 24.09784, 24.32606, 24.61293, 24.8998, 25.18667, 25.47354, 25.7604, 26.04727, 26.33414, 26.62101, 26.90788, 27.19475, 27.48162, 27.76848, 28.05535, 28.34222, 28.62909, 28.91596, 29.20283, 29.4897, 29.77657, 30.06343]) #[MeV]
Yield_86Sr_p_n_86Y_Uddin = np.array([0, 0.21684718, 1.652312932, 4.407782057, 8.760047767, 14.93860826, 23.12202572, 33.34681534, 45.62304025, 59.98849684, 76.35704094, 94.74345411, 115.0509768, 137.2021276, 161.0750919, 186.5016731, 213.2671419, 241.1139165, 269.7444951, 298.8292989, 328.0150147, 356.9347002, 385.2192722, 412.5097916, 438.4703389, 462.8008137, 485.2493291, 505.6232517, 523.7997569, 539.7318378, 553.4535178, 565.0802696, 574.8051198, 582.8894756, 589.6477613, 595.4250297, 600.566574, 605.3786123, 610.0790159, 614.5329708, 618.7000058, 621.5500839, 623.5937296, 625.7827022, 627.9007373, 629.9628781, 631.9827832, 633.9728768, 635.9443201, 637.9069907, 639.8696328, 641.8398722, 643.8242666, 645.8283223, 647.8566501, 649.912951, 651.9999972, 654.1197918, 656.2735891, 658.4619481, 660.6847883, 662.9414043, 665.2306358]) #[MBq/uA h]



path = './Yields_IAEA/'
filenames = ['sr6p86yt.txt', 'sr8p86yt.txt', 'rb5a86yt.txt']
# colors = ['mediumaquamarine', 'mediumorchid', 'deepskyblue']
# linestyles = [(0,(7,6)), ':', (0,(5,4,1,4))]
# # linestyles = ['--', ':', '-.']
# linewidths = [1, 1.7, 1.2]
linestyles = [(0,(10,6)), (0,(5,6)), (0,(5,4,1,4,1,4)), ':', 'solid']
linewidths = [1, 1, 1.3, 2, 1] 
colors = ['mediumorchid', 'mediumaquamarine', 'cornflowerblue', 'deepskyblue', 'grey']

i = 0
for filename in filenames:

    with open(path+filename, 'r') as file:
        lines = file.readlines()

        lines = [line.strip() for line in lines if line.strip()] #Remove empty lines
        label = lines[0] #First line gives the label for plotting
        data_lines = lines[4:] #choosing the lines with the yield data

        # Define lists for Energy og Yield
        energy_list = [0] #MeV
        yield_list = [0] #[MBq/uAh]

        # Prosesser hver data-linje
        for line in data_lines:
            columns = line.split()
            energy = float(columns[0])
            mbq_per_uah = float(columns[3])
            energy_list.append(energy)
            yield_list.append(mbq_per_uah)
    #1 MBq/C = 0.0036 MBq/uAh
    print(label)
    yields = np.array(yield_list)/0.0036 
    label_formatted = fr'$^{{{label[:3]}}}${label[3:-3]}$^{{{label[-3:-1]}}}${label[-1]}, IAEA'
    plt.plot(energy_list, yields, label = label_formatted, color = colors[i], linestyle = linestyles[i], linewidth=linewidths[i])
    i += 1

plt.plot(E_86Sr_p_n_86Y_Uddin, Yield_86Sr_p_n_86Y_Uddin/0.0036 , label = r'$^{86}$Sr(p,n)$^{86}$Y, Uddin (2025) [13]', color = colors[i], linestyle = linestyles[i], linewidth=linewidths[i])
i += 1
plt.plot(E_86Sr_d_2n_86Y_Uddin, Yield_86Sr_d_2n_86Y_Uddin/0.0036 , label = r'$^{86}$Sr(d,2n)$^{86}$Y, Uddin (2025) [13]', color = colors[i], linestyle = linestyles[i], linewidth=linewidths[i])

E_array_my_data =  np.array([0, 31.04155407, 36.351075, 41.11050957, 47.63516853])
# yield_list_my_data = np.array([0, 3.26210855, 17.62948948, 41.67509857, 80.41320086]) # (MBq/uA h)
yield_list_my_data = np.array([0,906.14126376,4897.08041243,11576.41626943,22337.00023931])


# plt.plot(energy_list, np.array(yield_list)*2, label = 'natZr(d,x)86Y', color = 'hotpink')
plt.plot(E_array_my_data, yield_list_my_data, label = r'$^{\text{nat}}$Zr(d,x)$^{86}$Y', color = 'hotpink', linewidth=2)
plt.xlabel('Beam energy (MeV)', fontsize=14)
# plt.ylabel(r'Physical Yields (MBq/$\mu$Ah)', fontsize=12)
# plt.ylabel('Physical Yields (MBq/C)', fontsize=14)
plt.ylabel('Physical Thick Target Yield (MBq/C)', fontsize=14)#,labelpad=7,fontsize=22)
plt.xlim(0,50)
plt.yscale('log')
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend(fontsize=11)
plt.savefig('./Figures/PhD/Physical_Yields_compare_routs_log.pdf',dpi=600)
plt.show()