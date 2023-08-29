
import numpy as np
import pandas as pd



data = pd.read_csv(r'CR060317_Zr04_18cm_30MeV_peak_data.csv')   
df = pd.DataFrame(data, columns=['isotope', 'energy'])


peakEnergy = 1462 #keV
peakEnergyUnc = 5 #keV

energies = df['energy']
isotopes = df['isotope']

for i in range(len(energies)):
    for j in range(len(energies)):
        sumEnergy = energies[i] + energies[j]

        # print(sumEnergy - peakEnergy)
        if (np.abs(sumEnergy - peakEnergy) < peakEnergyUnc):
            if (isotopes[i] == isotopes[j]):
                print (sumEnergy, energies[i], energies[j], isotopes[i], isotopes[j])