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
 
    sp.isotopes = ['43SC', '44SC', '44SCm', '46SC', '47SC', '47CA', '48V', '48SC']
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    # sp.plot(scale = 'linlin')
    sp.plot()
    sp.saveas(f'../../MyGeneratedFiles/Ti_foils/{spectrumName}/{spectrumName}_peak_data.csv')
   



if __name__ == '__main__':

    #Ti01:
    fit_peaks('CS060317_Ti01_18cm_30MeV', '18cm') #Date: 03/06/2017 10:23:02, $MEAS_TIM: 15386 15401

