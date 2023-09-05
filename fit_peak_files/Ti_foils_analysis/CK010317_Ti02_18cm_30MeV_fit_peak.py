import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):

    path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/30MeV/'
    spectrumFile = path + spectrumName + '.Spe'

    # cb = ci.Calibration('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_' + position + '_new.json')
    cb = ci.Calibration('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_' + position + '.json')
    # cb.plot()

    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 
 
    sp.isotopes = ['43SC', '44SC', '44SCm', '46SC', '47SC', '47CA', '48V', '48SC']
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    # sp.plot(scale = 'linlin')
    sp.plot()
    sp.saveas(f'/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/{spectrumName}/{spectrumName}_peak_data.csv')
   



if __name__ == '__main__':
  

    #Ti02:
    fit_peaks('CK010317_Ti02_18cm_30MeV', '18cm') #Date: 03/01/2017 15:43:23, $MEAS_TIM: 891 896

