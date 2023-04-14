import numpy as np 
import curie as ci 
import pandas as pd
import matplotlib.pyplot as plt


def calibration_18cm():

    cb = ci.Calibration()

    # sp_Ba133 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Ba133_150217_18cm.Spe')
    sp_Ba133 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Ba133_150217_18cm_mod.Spe')
    sp_Ba133.isotopes = ['133BA'] #['133BA', '56CO', '137CS', '152EU']
    sp_Ba133.auto_calibrate()
    # sp_Ba133.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/ElisesGeneratedFiles/DebugCalibration/Ba133.Spe')
    # sp_Ba133.plot()

    sp_Cs137 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Cs137_230317_18cm.Spe')
    sp_Cs137.isotopes = ['137CS'] 
    sp_Cs137.auto_calibrate()
    # sp_Cs137.plot()

    sp_newCs137 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/newCs137_030317_18cm.Spe')
    sp_newCs137.isotopes = ['137CS'] 
    sp_newCs137.auto_calibrate()
    # sp_newCs137.plot()

    sp_Cs137_1 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Cs137_230317_18cm(1).Spe')
    sp_Cs137_1.isotopes = ['137CS'] 
    sp_Cs137_1.auto_calibrate()
    # sp_Cs137_1.plot()

    sp_Cs137_2 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Cs137_240217_18cm.Spe')
    sp_Cs137_2.isotopes = ['137CS'] 
    sp_Cs137_2.auto_calibrate()
    # sp_Cs137_2.plot()

    sp_Eu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Eu152_150217_18cm.Spe')
    sp_Eu152.isotopes = ['152EU'] 
    sp_Eu152.auto_calibrate()
    # sp_Eu152.plot()

    sp_newEu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/newEu152_030317_18cm.Spe')
    sp_newEu152.isotopes = ['152EU'] 
    sp_newEu152.auto_calibrate()
    # sp_newEu152.plot()
    
    # sp_nuclearDataSchool = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/SommerskoleUSA2022/KevinsWork/IDM/calibration/AA20220307_EuCsBa_65cm_IDMV1.Spe')
    sp_nuclearDataSchool = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/SommerskoleUSA2022/KevinsWork/IDM/calibration/AA20220307_EuCsBa_65cm_IDMV1_mod.Spe')
    sp_nuclearDataSchool.isotopes = ['133BA'] 
    sp_nuclearDataSchool.auto_calibrate()
    sp_nuclearDataSchool.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/ElisesGeneratedFiles/DebugCalibration/nuclearDataSchool.Spe')
    # sp_nuclearDataSchool.plot()


    # sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
    #            {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
    #            {'isotope':'152EU', 'A0':3.929E4, 'ref_date':'01/01/2009 12:00:00'}]

    # sources = [{'isotope':'152EU', 'A0':3.929E4, 'ref_date':'01/01/2009 12:00:00'},
    #            {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'}]

    # sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
    #            {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'}]

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'}]

    # sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'}]

    sources = pd.DataFrame(sources)

    # cb.calibrate([sp_Ba133, sp_newCs137, sp_Cs137_1, sp_Cs137_2, sp_Eu152, sp_newEu152, sp_Cs137], sources=sources)
    # cb.calibrate([sp_Ba133, sp_Cs137, sp_newCs137, sp_Cs137_1], sources=sources)
    # cb.calibrate([sp_Cs137, sp_newCs137, sp_Cs137_1], sources=sources)
    # cb.calibrate([sp_Ba133], sources=sources)
    cb.calibrate([sp_Ba133], sources=sources)


    cb.plot()
    # cb.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/ElisesGeneratedFiles/Calibration/calibration_18cm.json')



def calibration_50cm():

    cb = ci.Calibration()

    sp_Ba133 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Ba133_090317_50cm.Spe')
    sp_Ba133.isotopes = ['133BA'] 
    sp_Ba133.auto_calibrate()
    sp_Ba133.plot()

    sp_Cs137 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Cs137_090317_50cm.Spe')
    sp_Cs137.isotopes = ['137CS'] 
    sp_Cs137.auto_calibrate()
    sp_Cs137.plot()

    sp_Eu152 = ci.Spectrum('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/calibration/Eu152_090317_50cm.Spe')
    sp_Eu152.isotopes = ['137CS'] 
    sp_Eu152.auto_calibrate()
    sp_Eu152.plot()

    sources = [{'isotope':'133BA', 'A0':3.989E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'137CS', 'A0':3.855E4, 'ref_date':'01/01/2009 12:00:00'},
               {'isotope':'152EU', 'A0':3.929E4, 'ref_date':'01/01/2009 12:00:00'}]

    sources = pd.DataFrame(sources)

    cb.calibrate([sp_Ba133, sp_Cs137, sp_Eu152], sources=sources)
    cb.plot()
    cb.saveas('/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/FilerFraAndrew/2017_Feb_Zr/ElisesGeneratedFiles/Calibration/calibration_50cm.json')




if __name__ == '__main__':
    calibration_18cm()
    # calibration_50cm()












