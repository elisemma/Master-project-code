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
 
    sp.isotopes = ['52MN', '54MN', '59FE', '55CO', '56CO', '57CO', '58CO', '58COm', '60CO','56NI', '57NI', '65NI', '60CU', '61CU', '64CU']
    # sp.isotopes = ['56CO', '58CO', '61CU']
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    sp.plot(scale = 'linlin')
    # sp.plot()
    sp.saveas(f'../../MyGeneratedFiles/Ni_foils/{spectrumName}/{spectrumName}_peak_data.csv')


if __name__ == '__main__':

    #Ni03:
    fit_peaks('DC200317_Ni03_18cm_30MeV', '18cm') #Date: 03/20/2017 15:39:25, $MEAS_TIM: 10045 10049


