import numpy as np 
import curie as ci 
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pl 




def calibration_10cm():
    cb = ci.Calibration()

    sp_Eu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Eu152_090317_10cm.Spe')
    sp_Eu152.isotopes = ['152EU']
    sp_Eu152.plot()

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'},
               {'isotope':'56CO', 'A0':3.929E4, 'ref_date':'01/01/2009 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Eu152], sources=sources)
    cb.plot(show=True, saveas = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/Figures/calibration_plots_10cm.pdf')
    cb.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_10cm.json')



def calibration_18cm():
    cb = ci.Calibration()

    sp_Ba133 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Ba133_150217_18cm.Spe')
    sp_Ba133.isotopes = ['133BA'] 
    # sp_Ba133.plot()

    sp_Co56 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Co56_190217_18cm.Spe')
    sp_Co56.isotopes = ['56CO'] 
    # sp_Co56.plot()

    sp_Cs137 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Cs137_230317_18cm.Spe')
    sp_Cs137.isotopes = ['137CS'] 
    # sp_Cs137.plot()

    sp_Cs137_2 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Cs137_240217_18cm.Spe')
    sp_Cs137_2.isotopes = ['137CS'] 
    # sp_Cs137_2.plot()

    sp_Eu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Eu152_150217_18cm.Spe')
    sp_Eu152.isotopes = ['152EU'] 
    # fig = plt.figure()
    sp_Eu152.plot()  
    # pl.dump(fig,open('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/152Eu18cm.pickle','wb'))  

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'},
               {'isotope':'56CO', 'A0':3.929E4, 'ref_date':'01/01/2009 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Ba133, sp_Cs137, sp_Cs137_2, sp_Eu152, sp_Co56], sources=sources)
    cb.plot()
    # cb.plot(show=False, saveas = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/Figures/calibration_plots_18cm.pdf')
    # cb.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_18cm.json')



def calibration_18cm_new():
    cb = ci.Calibration()

    sp_newCs137 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/newCs137_030317_18cm.Spe')
    sp_newCs137.isotopes = ['137CS'] 
    # sp_newCs137.plot()

    sp_newEu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/newEu152_030317_18cm.Spe')
    sp_newEu152.isotopes = ['152EU'] 
    # sp_newEu152.plot()

    sources = [{'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_newCs137, sp_newEu152], sources=sources)
    cb.plot(show=False, saveas = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/Figures/calibration_plots_18cm_new.pdf')
    cb.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_18cm_new.json')
    # cb.plot()


def calibration_40cm():
    cb = ci.Calibration()

    sp_Eu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Eu152_150217_40cm.Spe')
    sp_Eu152.isotopes = ['152EU']
    # sp_Eu152.plot()

    sources = [{'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Eu152], sources=sources)
    cb.plot(show=False, saveas = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/Figures/calibration_plots_40cm.pdf')
    cb.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_40cm.json')



def calibration_50cm():

    cb = ci.Calibration()

    sp_Ba133 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Ba133_090317_50cm.Spe')
    sp_Ba133.isotopes = ['133BA'] 
    # sp_Ba133.plot()

    sp_Cs137 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Cs137_090317_50cm.Spe')
    sp_Cs137.isotopes = ['137CS'] 
    # sp_Cs137.plot()

    sp_Eu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Eu152_090317_50cm.Spe')
    sp_Eu152.isotopes = ['152EU'] 
    # sp_Eu152.plot()

    sp_Co56 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/2017_Feb_Zr/calibration/Co56_140317_50cm.Spe')
    sp_Co56.isotopes = ['56CO'] 
    # sp_Co56.plot()

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'152EU', 'A0':370000, 'ref_date':'11/01/1984 12:00:00'},
               {'isotope':'56CO', 'A0':3.929E4, 'ref_date':'01/01/2009 12:00:00'}]
    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Ba133, sp_Cs137, sp_Eu152, sp_Co56], sources=sources)
    cb.plot(show=False, saveas = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/Figures/calibration_plots_50cm.pdf')
    cb.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Calibration/json_files/calibration_50cm.json')





if __name__ == '__main__':

    calibration_10cm()
    # calibration_18cm()
    # calibration_18cm_new()
    # calibration_40cm()
    # calibration_50cm()













