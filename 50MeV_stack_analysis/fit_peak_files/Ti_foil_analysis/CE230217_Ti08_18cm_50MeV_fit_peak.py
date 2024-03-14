import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):
    
    path = '../../../2017_Feb_Zr/50MeV/'

    spectrumFile = path + spectrumName + '.Spe'

    cb = ci.Calibration('../../MyGeneratedFiles/Calibration/json_files/calibration_' + position + '.json')

    
    # cb.plot()

    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 
 
    sp.isotopes = ['43SC', '44SC', '44SCm', '46SC', '47SC', '47CA', '48V', '48SC']
    # sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0, 'skew_fit':True}
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}

    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    # sp.plot(scale = 'linlin')
    # sp.rebin(4000)
    sp.plot()
    
    sp.saveas(f'../../MyGeneratedFiles/Ti_foils/{spectrumName}/{spectrumName}_peak_data.csv')

   



if __name__ == '__main__':

    fit_peaks('CE230217_Ti08_18cm_50MeV', '18cm')