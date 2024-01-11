import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):

    path = '../../2017_Feb_Zr/30MeV/'
    spectrumFile = path + spectrumName + '.Spe'

    cb = ci.Calibration('../../MyGeneratedFiles/Calibration/json_files/calibrationByHand_DF240317_Ni04_18cm_30MeV.json')
    
    # cb.plot()

    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 
 
    # sp.isotopes = ['56CO', '58CO']
    sp.isotopes = ['52MN', '54MN', '59FE', '55CO', '56CO', '57CO', '58CO', '58COm', '60CO','56NI', '57NI', '65NI', '60CU', '61CU', '64CU']
    
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0, 'skew_fit':True}
    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    sp.plot(scale = 'linlin')
    #sp.plot()
    sp.saveas(f'../../MyGeneratedFiles/Ni_foils/{spectrumName}/{spectrumName}_peak_data.csv')
   



if __name__ == '__main__':

    #Ni04:
    fit_peaks('DF240317_Ni04_18cm_30MeV', '18cm') #Date: 03/23/2017 15:58:09, $MEAS_TIM: 68418 68435
