import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):

    path = '../../2017_Feb_Zr/30MeV/'
    spectrumFile = path + spectrumName + '.Spe'

    cb = ci.Calibration('../../MyGeneratedFiles/Calibration/json_files/calibration_' + position + '_new.json')
    
    # cb.plot()

    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 
 
    sp.isotopes = ['43K', '43SC', '44SC', '44SCm', '46SC', '47SC', '47CA', '48V', '48SC']
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    # sp.plot(scale = 'linlin')
    sp.plot()
    sp.saveas(f'../../MyGeneratedFiles/Ti_foils/{spectrumName}/{spectrumName}_peak_data.csv')
   



if __name__ == '__main__':

    #Ti05:
    fit_peaks('CP030317_Ti05_18cm_30MeV', '18cm') #Date: 03/02/2017 14:50:40, $MEAS_TIM: 72299 72303

