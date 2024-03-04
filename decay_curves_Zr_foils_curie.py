
import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import numpy as np 



# Zr peak data _____________________________________________________________________________________________________________________
path_Zr = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Zr_foils/'
file_AXZR05 = 'AX130217_Zr05_18cm_30MeV/AX130217_Zr05_18cm_30MeV_peak_data.csv'
file_BAZR01 = 'BA130217_Zr01_18cm_30MeV/BA130217_Zr01_18cm_30MeV_peak_data.csv'
file_BIZR04 = 'BI140217_Zr04_18cm_30MeV/BI140217_Zr04_18cm_30MeV_peak_data.csv'
file_BPZR01 = 'BP150217_Zr01_18cm_30MeV/BP150217_Zr01_18cm_30MeV_peak_data.csv'
file_BRZR02 = 'BR150217_Zr02_18cm_30MeV/BR150217_Zr02_18cm_30MeV_peak_data.csv'
file_BSZR05 = 'BS160217_Zr05_18cm_30MeV/BS160217_Zr05_18cm_30MeV_peak_data.csv'
file_BTZR03 = 'CH260217_Zr03_18cm_30MeV/CH260217_Zr03_18cm_30MeV_peak_data.csv'
file_BXZR01 = 'BX190217_Zr01_18cm_30MeV/BX190217_Zr01_18cm_30MeV_peak_data.csv'
file_CHZR03 = 'BT160217_Zr03_18cm_30MeV/BT160217_Zr03_18cm_30MeV_peak_data.csv'
file_CRZR04 = 'CR060317_Zr04_18cm_30MeV/CR060317_Zr04_18cm_30MeV_peak_data.csv'

file_BFZR01 = 'BF130217_Zr01_40cm_30MeV/BF130217_Zr01_40cm_30MeV_peak_data.csv'
file_BGZR02 = 'BG130217_Zr02_40cm_30MeV/BG130217_Zr02_40cm_30MeV_peak_data.csv'
file_BHZR03 = 'BH130217_Zr03_40cm_30MeV/BH130217_Zr03_40cm_30MeV_peak_data.csv'
file_BOZR04 = 'BO150217_Zr04_40cm_30MeV/BO150217_Zr04_40cm_30MeV_peak_data.csv'

file_CYZR02 = 'CY090317_Zr02_10cm_30MeV/CY090317_Zr02_10cm_30MeV_peak_data.csv'
file_CZZR05 = 'CZ120317_Zr05_10cm_30MeV/CZ120317_Zr05_10cm_30MeV_peak_data.csv'

file_concat_Zr = 'combined_peak_data_Zr.csv'

df_AXXR05 = pd.read_csv(path_Zr+file_AXZR05)
df_BAZR01 = pd.read_csv(path_Zr+file_BAZR01)
df_BIZR04 = pd.read_csv(path_Zr+file_BIZR04)
df_BPZR01 = pd.read_csv(path_Zr+file_BPZR01)
df_BRZR02 = pd.read_csv(path_Zr+file_BRZR02)
df_BSZR05 = pd.read_csv(path_Zr+file_BSZR05)
df_BTZR03 = pd.read_csv(path_Zr+file_BTZR03)
df_BXZR01 = pd.read_csv(path_Zr+file_BXZR01)
df_CHZR03 = pd.read_csv(path_Zr+file_CHZR03)
df_CRZR04 = pd.read_csv(path_Zr+file_CRZR04)

df_BFZR01 = pd.read_csv(path_Zr+file_BFZR01)
df_BGZR02 = pd.read_csv(path_Zr+file_BGZR02)
df_BHZR03 = pd.read_csv(path_Zr+file_BHZR03)
df_BOZR04 = pd.read_csv(path_Zr+file_BOZR04)

df_CYZR02 = pd.read_csv(path_Zr+file_CYZR02)
df_CZZR05 = pd.read_csv(path_Zr+file_CZZR05)



df_concat_Zr = pd.concat((df_AXXR05, df_BAZR01, df_BIZR04, df_BPZR01, df_BRZR02, df_BSZR05, df_BTZR03, df_BXZR01, df_CHZR03, df_CRZR04, df_BFZR01, df_BGZR02, df_BHZR03, df_BOZR04, df_CYZR02, df_CZZR05), axis = 0)
df_concat_Zr= df_concat_Zr[(df_concat_Zr['isotope'] != '90NB') | (df_concat_Zr['energy'] != 329.058)] # Exclude rows where isotope is '90NB' and energy is 329.058
df_concat_Zr.to_csv(path_Zr+file_concat_Zr)





def fit_prod_rate(isotope_list, isotope_chain_parent, foil, path, file, plot=False):
    t_irr_h = 0.33

    R = {}
    for isotope in isotope_list:
       R[isotope] =[[1E4,t_irr_h]] 

    dc = ci.DecayChain(isotope_chain_parent, R=R, units='h')
    dc.get_counts(foil, '02/13/2017 14:47:00', path+file)
    isotopes, R, cov_R = dc.fit_R()
    print(f'A0 for {isotope_chain_parent}: ', dc.activity(isotope_chain_parent, 0))

    if plot==True:
        dc.plot()

    return isotopes, R, cov_R




def calc_prod_rates_in_foil(isotope_list_list, isotope_chain_parent_list, foil, write_to_file=False, show_plot=False):
    print(' ')
    print(foil)
    csv_file_path = f'./Calculated_R/{foil}_R_by_curie.csv'
    if write_to_file==True:
            with open(csv_file_path, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerow(['Isotope', 'R', 'R_unc', 'Rel_R_unc_in_percent'])

    for isotope_list, parent in zip(isotope_list_list, isotope_chain_parent_list):
        isotopes, R, cov_R = fit_prod_rate(isotope_list, parent, foil, path_Zr, file_concat_Zr, plot=show_plot)

        print(isotopes)
        print(R)
        # print(cov_R)
        if write_to_file==True:
            with open(csv_file_path, 'a', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                for i, isotope in enumerate(isotope_list):
                    # print(i, isotope)
                    csv_writer.writerow([f'{isotope}', f'{R[i]}', f'{np.sqrt(cov_R[i,i])}', f'{np.sqrt(cov_R[i,i])/R[i]*100:.2e}'])




       
# isotope_list_list_Zr01 = [['87NB', '87Ym', '87Y'], ['89NBm', '89ZR'], ['90NB'], ['92Y'], ['95NBm', '95NB'], ['96NB']]
# isotope_chain_parent_list_Zr01 = ['87NB', '89NBm', '90NB', '92Y', '95NBm', '96NB']

# isotope_list_list_Zr02 = [['87NB', '87Ym', '87Y'], ['89NBm', '89ZR'], ['90NB', '90Ym'], ['92Y'], ['95NBm', '95NB', '95ZR'], ['96NB']]
# isotope_chain_parent_list_Zr02 = ['87NB', '89NBm', '90NB', '92Y', '95NBm', '96NB']

# isotope_list_list_Zr03 = [['87NB'], ['89NBm', '89ZR'], ['90NB'], ['91Y'], ['92Y'], ['95NBm', '95NB', '95ZR'], ['96NB']]
# isotope_chain_parent_list_Zr03 = ['87NB', '89NBm', '90NB', '91Y', '92Y', '95NBm', '96NB']

# isotope_list_list_Zr04 = [['88Y'], ['89NB', '89Ym', '89ZR'], ['90NB'], ['91Y'], ['92Y'], ['95NBm', '95NB', '95ZR'], ['96NB']]
# isotope_chain_parent_list_Zr04 = ['88Y', '89NB', '90NB', '91Y', '92Y', '95NBm', '96NB']

# isotope_list_list_Zr05 = [['89NB', '89ZR'], ['90NB', '90Ym'], ['92Y'], ['95NB', '95ZR']]
# isotope_chain_parent_list_Zr05 = ['89NB', '90NB', '92Y', '95NB']


isotope_list_list_Zr01 = [['90NB'], ['96NB']]
isotope_chain_parent_list_Zr01 = ['90NB', '96NB']

isotope_list_list_Zr02 = [['90NB'], ['96NB']]
isotope_chain_parent_list_Zr02 = ['90NB', '96NB']

isotope_list_list_Zr03 = [['90NB'], ['96NB']]
isotope_chain_parent_list_Zr03 = ['90NB', '96NB']

isotope_list_list_Zr04 = [['90NB'], ['96NB']]
isotope_chain_parent_list_Zr04 = ['90NB', '96NB']

isotope_list_list_Zr05 = [['90NB'], ['96NB']]
isotope_chain_parent_list_Zr05 = ['90NB', '96NB']


# isotope_list_list_Zr01 = [['87NB', '87Ym', '87Y']]
# isotope_chain_parent_list_Zr01 = ['87NB']

# isotope_list_list_Zr02 = [['87NB', '87Ym', '87Y']]
# isotope_chain_parent_list_Zr02 = ['87NB']

# # isotope_list_list_Zr03 = [['87NB']]
# # isotope_chain_parent_list_Zr03 = ['87NB']

# isotope_list_list_Zr04 = [['88Y'], ['89NB', '89Ym', '89ZR'], ['90NB'], ['91Y'], ['92Y'], ['95NBm', '95NB', '95ZR'], ['96NB']]
# isotope_chain_parent_list_Zr04 = ['88Y', '89NB', '90NB', '91Y', '92Y', '95NBm', '96NB']

# isotope_list_list_Zr05 = [['89NB', '89ZR'], ['90NB', '90Ym'], ['92Y'], ['95NB', '95ZR']]
# isotope_chain_parent_list_Zr05 = ['89NB', '90NB', '92Y', '95NB']




calc_prod_rates_in_foil(isotope_list_list_Zr01, isotope_chain_parent_list_Zr01, 'Zr01', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Zr02, isotope_chain_parent_list_Zr02, 'Zr02', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Zr03, isotope_chain_parent_list_Zr03, 'Zr03', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Zr04, isotope_chain_parent_list_Zr04, 'Zr04', write_to_file=True, show_plot=False)
calc_prod_rates_in_foil(isotope_list_list_Zr05, isotope_chain_parent_list_Zr05, 'Zr05', write_to_file=True, show_plot=False)














































# t_irr_h = 0.33
# R = {'87NB':[[1E4,1]],
#          '87Y':[[1E2,1]],
#          '87Ym':[[1E4,1]]}
# dc = ci.DecayChain('87NB', R=R, units='h')
# dc.get_counts('Zr01', '02/13/2017 14:27:00', path_Zr+file_concat_Zr)
# isotopes, R, cov_R = dc.fit_R()


# # dc = ci.DecayChain('87NB', R={'87NB':[[3E5,1.0]], '87Y':[[1E6,1.0]], '87Ym':[[1E6,1.0]]}, units='h')
# # dc.get_counts('Zr01', '02/13/2017 14:27:00', 'peak_data.csv')

# # isotopes, R, var_R = dc.fit_R()
# # dc.plot(palette='american', shade='dark', style='paper', max_plot_chi2=100)

# print(isotopes)
# print(R)
# print(cov_R)
# print(f'{np.sqrt(cov_R[0,0]):.3e}')
# # print(cov_R)

# # A0_Nb87 = dc.activity('87NB', 0)

# # A0_Y87 = dc.activity('87Y', 0)
# # A0_Y87m = dc.activity('87Ym',0)

# # print(A0_Nb87, A0_unc_Nb87)
# # print(A0_unc_Nb87/A0_Nb87*100)
# # print(A0_Y87)
# # print(A0_Y87m)
