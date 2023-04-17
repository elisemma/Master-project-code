import numpy as np 
import pandas as pd
import curie as ci
import matplotlib.pyplot as plt 



def fit_peaks(spectrumName, position):

    path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/30MeV/'
    spectrumFile = path + spectrumName + '.Spe'

    cb = ci.Calibration('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/ElisesGeneratedFiles/Calibration/calibration_' + position + '.json')
    # cb.plot()

    sp = ci.Spectrum(spectrumFile)
    sp.cb = cb 

    effcal = list(cb.effcal)
    cb.effcal = effcal

    sp.isotopes = ['86Y', '87Y', '88Y', '88ZR', '87NB', '88NB', '89NB', '90NB', '87MO', '89MO']
    sp.fit_peaks(bg='constant')
    sp._peaks = sp.peaks[(sp.peaks['counts']>300)]
    sp._peaks = sp.peaks[(sp.peaks['chi2']<10)]
    sp.plot()





if __name__ == '__main__':

    fit_peaks('BA130217_Zr01_18cm_30MeV', '18cm')
