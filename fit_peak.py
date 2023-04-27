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

    # sp.isotopes = ['89NB', '90NB', '91NB', '92NB', '93NB', '94NB', '95NB', '96NB', '97NB', '98NB', 
    #                '90ZR', '91ZR', '92ZR', '93ZR', '94ZR', '95ZR', '96ZR', 
    #                '86Y', '87Y', '88Y', '89Y', '90Y', '91Y', '92Y', '93Y', '94Y']

    # sp.isotopes = ['89NB']

    # peaks = {'energy': [511, 769.6], 'intensity': [162, 6.2], 'unc_intensity': [17, 0.6]}
    # peaks = pd.DataFrame(peaks)

    # peaks = {'energy': 1611, 'intensity': 100, 'unc_intensity': 10}
    # # peaks = pd.DataFrame(peaks)


    # # sp.fit_peaks(gammas = [peak1, peak2], bg='constant')
    # # sp.fit_config = {'bg':'constant'}
    # # sp.fit_config = {'bg': 'linear', 'gammas': peaks}
    # sp.fit_peaks(gammas=[{'energy':1611, 'intensity':100, 'unc_intensity':10}])

    # sp._peaks = sp.peaks[(sp.peaks['counts']>300)]
    # sp._peaks = sp.peaks[(sp.peaks['chi2']<10)]
    sp.isotopes = ['89NB']
    sp.fit_peaks(gammas=[{'energy':1100.0, 'intensity':10.66, 'unc_intensity':0.55}]) #E=1100 works, E=1807 does not work
    # sp.fit_peaks(gammas=ci.Isotope('40K').gammas(istp_col=True))
    sp.plot()


    """
    Eksempel:
    sp.isotopes = ['152EU', '40K']
    sp.fit_peaks(gammas=[{'energy':1460.8, 'intensity':10.66, 'unc_intensity':0.55}])
    sp.fit_peaks(gammas=ci.Isotope('40K').gammas(istp_col=True))
    sp.plot()
    """




if __name__ == '__main__':

    fit_peaks('BA130217_Zr01_18cm_30MeV', '18cm')
