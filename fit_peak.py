import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):

    path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/30MeV/'
    spectrumFile = path + spectrumName + '.Spe'

    cb = ci.Calibration('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_' + position + '.json')
    # cb.plot()

    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 
    # sp.isotopes = ['87NB', '89NB', '90NB', '91NB', '92NB', '94NB', '95NB', '96NB', '97NB', '98NB', 
    #                '90ZR', '95ZR',
    #                '86Y', '87Y', '92Y', '94Y',

    #                '26AL',
    #                '89NBm', '94NBm', '97NBm', '98NBm', 
    #                '86Ym', '87Ym', '89Ym', '90Ym', '91Ym']#93Nb, 91,92,93,94,96zr, 88,89,90,91,93Y har for lang half life/ er stabil, 90Y har ikke såå lang half life, men man ser den kanskje senere, 94Y burde vi vel ha sett?, hvorfor ser jeg ikke 86mY og 91mY


    sp.isotopes = ['86Y', '87NB', '87Y', '87Ym', '89NB', '89NBm', '89Ym', '89ZR', '90NB', '90Ym', '92Y', '95NB', '95NBm', '96NB', '97NB', '97NBm', '98NBm']
    # sp.isotopes = ['97ZR']

    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}
    sp.plot(scale = 'linlin')
    # sp.plot()
   



if __name__ == '__main__':

    #Zr01:
    # fit_peaks('BA130217_Zr01_18cm_30MeV', '18cm') #Date: 02/13/2017 20:36
    # fit_peaks('BP150217_Zr01_18cm_30MeV', '18cm') #Date: 02/15/2017 10:39
    fit_peaks('BX190217_Zr01_18cm_30MeV', '18cm') #Date: 02/17/2017 18:18
    # #Zr02:
    # fit_peaks('BR150217_Zr02_18cm_30MeV', '18cm') #Date: 02/15/2017 19:13
    # #Zr03:
    # fit_peaks('BT160217_Zr03_18cm_30MeV', '18cm') #Date: 02/16/2017 10:32
    # fit_peaks('CH260217_Zr03_18cm_30MeV', '18cm')
    # #Zr04:
    # fit_peaks('BI140217_Zr04_18cm_30MeV', '18cm')
    # fit_peaks('CR060317_Zr04_18cm_30MeV', '18cm')
    # #Zr05:
    # fit_peaks('AX130217_Zr05_18cm_30MeV', '18cm')
    # fit_peaks('BS160217_Zr05_18cm_30MeV', '18cm')





