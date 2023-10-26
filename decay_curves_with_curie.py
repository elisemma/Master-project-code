import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def fit_prod_rate(isotope, foil, path, file):
    R_estimated = 5
    t_irr_h = 0.33
    dc = ci.DecayChain(isotope, units='h', R=[[R_estimated, t_irr_h]])
    dc.get_counts(foil, '02/13/2017 14:27:00', path+file)

    isotopes, R, var_R = dc.fit_R()
    dc.plot(titel = foil)



path_Ti = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Ti_foils/'
file_CJTi01 = 'CJ010317_Ti01_18cm_30MeV/CJ010317_Ti01_18cm_30MeV_peak_data.csv'
file_CKTi02 = 'CK010317_Ti02_18cm_30MeV/CK010317_Ti02_18cm_30MeV_peak_data.csv'
file_CLTi03 = 'CL010317_Ti03_18cm_30MeV/CL010317_Ti03_18cm_30MeV_peak_data.csv'
file_CMTi04 = 'CM010317_Ti04_18cm_30MeV/CM010317_Ti04_18cm_30MeV_peak_data.csv'
file_CPTi05 = 'CP030317_Ti05_18cm_30MeV/CP030317_Ti05_18cm_30MeV_peak_data.csv'
file_CQTi04 = 'CQ030317_Ti04_18cm_30MeV/CQ030317_Ti04_18cm_30MeV_peak_data.csv'
file_CSTi01 = 'CS060317_Ti01_18cm_30MeV/CS060317_Ti01_18cm_30MeV_peak_data.csv'
file_CTTi02 = 'CT060317_Ti02_18cm_30MeV/CT060317_Ti02_18cm_30MeV_peak_data.csv'
file_CUTi03 = 'CU060317_Ti03_18cm_30MeV/CU060317_Ti03_18cm_30MeV_peak_data.csv'

file_concat_Ti = 'combined_peak_data_Ti.csv'

df_CJTi01 = pd.read_csv(path_Ti+file_CJTi01)
df_CKTi02 = pd.read_csv(path_Ti+file_CKTi02)
df_CLTi03 = pd.read_csv(path_Ti+file_CLTi03)
df_CMTi04 = pd.read_csv(path_Ti+file_CMTi04)
df_CPTi05 = pd.read_csv(path_Ti+file_CPTi05)
df_CQTi04 = pd.read_csv(path_Ti+file_CQTi04)
df_CSTi01 = pd.read_csv(path_Ti+file_CSTi01)
df_CTTi02 = pd.read_csv(path_Ti+file_CTTi02)
df_CUTi03 = pd.read_csv(path_Ti+file_CUTi03)

df_concat = pd.concat((df_CJTi01, df_CKTi02, df_CLTi03, df_CMTi04, df_CPTi05, df_CQTi04, df_CSTi01, df_CTTi02, df_CUTi03), axis = 0)
df_concat.to_csv(path_Ti+file_concat_Ti)

isotope = '48V'
foil = 'Ti05'
path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/'
file = 'combined_peak_data_Ti01.csv'

fit_prod_rate(isotope, foil, path_Ti, file_concat_Ti)











