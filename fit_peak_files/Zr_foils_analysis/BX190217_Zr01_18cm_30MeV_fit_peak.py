import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):

    path = '../../2017_Feb_Zr/30MeV/'
    spectrumFile = path + spectrumName + '.Spe'

    cb = ci.Calibration('../../MyGeneratedFiles/Calibration/json_files/calibration_' + position + '.json')
    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 
   
    sp.isotopes = ['87Y', '87Ym', '87SR', '89NB', '89ZR', '90NB', '90Ym', '91NBm',  '92NBm', '95NB', '95NBm', '24NA', '40K', '96NB', '88Y']
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    # sp.plot(scale = 'linlin')
    sp.plot()
    sp.saveas(f'../../MyGeneratedFiles/Zr_foils/{spectrumName}/{spectrumName}_peak_data.csv')
   



if __name__ == '__main__':

    #Zr01:
    # fit_peaks('BA130217_Zr01_18cm_30MeV', '18cm') #Date: 02/13/2017 20:36
    # fit_peaks('BP150217_Zr01_18cm_30MeV', '18cm') #Date: 02/15/2017 10:39
    fit_peaks('BX190217_Zr01_18cm_30MeV', '18cm') #Date: 02/17/2017 18:18
    # #Zr02:
    # fit_peaks('BR150217_Zr02_18cm_30MeV', '18cm') #Date: 02/15/2017 19:13
    # #Zr03:
    # fit_peaks('BT160217_Zr03_18cm_30MeV', '18cm') #Date: 02/16/2017 10:32
    # fit_peaks('CH260217_Zr03_18cm_30MeV', '18cm') #Date: 02/24/2017 17:15
    # #Zr04:
    # fit_peaks('BI140217_Zr04_18cm_30MeV', '18cm') #Date: 02/14/2017 08:49
    # fit_peaks('CR060317_Zr04_18cm_30MeV', '18cm') #Date: 03/03/2017 19:04  NB: Must use new calibration 
    # #Zr05:
    # fit_peaks('AX130217_Zr05_18cm_30MeV', '18cm') #Date: 02/13/2017 15:07
    # fit_peaks('BS160217_Zr05_18cm_30MeV', '18cm') #Date: 02/15/2017 21:24





