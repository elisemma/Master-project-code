
import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):

    path = '../../2017_Feb_Zr/30MeV/'
    spectrumFile = path + spectrumName + '.Spe'

    cb = ci.Calibration('../../MyGeneratedFiles/Calibration/json_files/calibration_' + position + '.json')
  
    # cb.plot()

    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 

    sp.isotopes = ['87NB', '87Y', '87Ym', '89NB', '89NBm', '89ZR', '90NB', '90Ym', '90Y', '92Y', '96NB', '24NA', '40K', '96NB']
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    # sp.plot(scale = 'linlin')
    sp.plot()
    sp.saveas(f'../../MyGeneratedFiles/Zr_foils/{spectrumName}/{spectrumName}_peak_data.csv')
   



if __name__ == '__main__':

  
    fit_peaks('BG130217_Zr02_40cm_30MeV', '40cm') #Date: 02/13/2017 15:07

