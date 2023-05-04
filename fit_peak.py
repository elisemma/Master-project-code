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
    sp.isotopes = ['87NB', '89NB', '90NB', '91NB', '92NB', '94NB', '95NB', '96NB', '97NB', '98NB', 
                   '90ZR', '95ZR',
                   '86Y', '87Y', '92Y', '94Y',

                   '26AL',
                   '89NBm', '94NBm', '97NBm', '98NBm', 
                   '86Ym', '87Ym', '89Ym', '90Ym', '91Ym']#93Nb, 91,92,93,94,96zr, 88,89,90,91,93Y har for lang half life/ er stabil, 90Y har ikke såå lang half life, men man ser den kanskje senere, 94Y burde vi vel ha sett?, hvorfor ser jeg ikke 86mY og 91mY
    # sp.fit_config = {'skew_fit': True}
    sp.plot()
   



if __name__ == '__main__':

    fit_peaks('BA130217_Zr01_18cm_30MeV', '18cm')
