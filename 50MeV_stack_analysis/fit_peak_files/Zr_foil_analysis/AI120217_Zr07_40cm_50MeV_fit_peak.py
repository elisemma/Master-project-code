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
 
    sp.isotopes = ['86Y', '86Ym', '87Y', '87Ym', '88Y', '89Y', '89Ym', '90Y', '90Ym', '91Y', '92Y',
                   '86NB', '87NB', '88NB', '89NB', '89NBm', '90NB', '95NB', '95NBm', '96NB',
                   '89ZR',
                   '24NA', '40K']
    # sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0, 'skew_fit':True}
    sp.fit_config = {'SNR_min':3.5, 'dE_511':9.0}

    # sp.fit_config = {'SNR_min':2, 'dE_511':9.0}
    # sp.plot(scale = 'linlin')
    # sp.rebin(4000)
    sp.plot()
    
    sp.saveas(f'../../MyGeneratedFiles/Zr_foils/{spectrumName}/{spectrumName}_peak_data.csv')

   



if __name__ == '__main__':

    fit_peaks('AI120217_Zr7_40cm_50MeV', '40cm')