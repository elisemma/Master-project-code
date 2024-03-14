
import curie as ci
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import numpy as np 



# Zr peak data _____________________________________________________________________________________________________________________
path_Zr = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Zr_foils/'
# 30 MeV:
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

# 50 MeV:
file_AAZR01 = 'AA120217_Zr10_50cm_50MeV/AA120217_Zr10_50cm_50MeV_peak_data.csv'
file_ABZR09 = 'AB120217_Zr09_50cm_50MeV/AB120217_Zr09_50cm_50MeV_peak_data.csv'
file_ABZR08 = 'AC120217_Zr08_50cm_50MeV/AC120217_Zr08_50cm_50MeV_peak_data.csv'
file_ACZR07 = 'AD120217_Zr07_50cm_50MeV/AD120217_Zr07_50cm_50MeV_peak_data.csv'
file_AEZR06 = 'AE120217_Zr06_50cm_50MeV/AE120217_Zr06_50cm_50MeV_peak_data.csv'
file_AFZR10 = 'AF120217_Zr10_40cm_50MeV/AF120217_Zr10_40cm_50MeV_peak_data.csv'
file_AGZR09 = 'AG120217_Zr9_40cm_50MeV/AG120217_Zr9_40cm_50MeV_peak_data.csv'
file_AHZR08 = 'AH120217_Zr8_40cm_50MeV/AH120217_Zr8_40cm_50MeV_peak_data.csv'
file_AIZR07 = 'AI120217_Zr7_40cm_50MeV/AI120217_Zr7_40cm_50MeV_peak_data.csv'
file_AJZR06 = 'AJ120217_Zr6_40cm_50MeV/AJ120217_Zr6_40cm_50MeV_peak_data.csv'
file_AKZR06 = 'AK120217_Zr6_15cm_50MeV/AK120217_Zr6_15cm_50MeV_peak_data.csv'
file_ALZR06 = 'AL130217_Zr6_18cm_50MeV/AL130217_Zr6_18cm_50MeV_peak_data.csv'
file_AMZR06 = 'AM130217_Zr6_22cm_50MeV/AM130217_Zr6_22cm_50MeV_peak_data.csv'
file_ANZR07 = 'AN130217_Zr7_22cm_50MeV/AN130217_Zr7_22cm_50MeV_peak_data.csv'
file_AOZR07 = 'AO130217_Zr7_40cm_50MeV/AO130217_Zr7_40cm_50MeV_peak_data.csv'
file_APZR08 = 'AP130217_Zr8_22cm_50MeV/AP130217_Zr8_22cm_50MeV_peak_data.csv'
file_AQZR09 = 'AQ130217_Zr9_22cm_50MeV/AQ130217_Zr9_22cm_50MeV_peak_data.csv'
file_ARZR10 = 'AR130217_Zr10_22cm_50MeV/AR130217_Zr10_22cm_50MeV_peak_data.csv'
file_ASZR07 = 'AS130217_Zr07_18cm_50MeV/AS130217_Zr07_18cm_50MeV_peak_data.csv'
file_AUZR08 = 'AU130217_Zr08_18cm_50MeV/AU130217_Zr08_18cm_50MeV_peak_data.csv'
file_AVZR09 = 'AV130217_Zr09_18cm_50MeV/AV130217_Zr09_18cm_50MeV_peak_data.csv'
file_AWZR10 = 'AW130217_Zr10_18cm_50MeV/AW130217_Zr10_18cm_50MeV_peak_data.csv'
file_BJZR06 = 'BJ140217_Zr06_10cm_50MeV/BJ140217_Zr06_10cm_50MeV_peak_data.csv'
file_BKZR07 = 'BK140217_Zr07_10cm_50MeV/BK140217_Zr07_10cm_50MeV_peak_data.csv'
file_BLZR08 = 'BL140217_Zr08_10cm_50MeV/BL140217_Zr08_10cm_50MeV_peak_data.csv'
file_BMZR09 = 'BM140217_Zr09_10cm_50MeV/BM140217_Zr09_10cm_50MeV_peak_data.csv'
file_BNZR10 = 'BN140217_Zr10_10cm_50MeV/BN140217_Zr10_10cm_50MeV_peak_data.csv'
file_BQZR06 = 'BQ150217_Zr06_18cm_50MeV/BQ150217_Zr06_18cm_50MeV_peak_data.csv'
file_CAZR06 = 'CA210217_Zr06_18cm_50MeV/CA210217_Zr06_18cm_50MeV_peak_data.csv'
file_CBZR07 = 'CB220217_Zr07_18cm_50MeV/CB220217_Zr07_18cm_50MeV_peak_data.csv'
file_CDZR08 = 'CD230217_Zr08_18cm_50MeV/CD230217_Zr08_18cm_50MeV_peak_data.csv'
file_CNZR09 = 'CN020317_Zr09_18cm_50MeV/CN020317_Zr09_18cm_50MeV_peak_data.csv'
file_CVZR10 = 'CV020317_Zr10_18cm_50MeV/CV020317_Zr10_18cm_50MeV_peak_data.csv'
file_CVZR10 = 'CV080317_Zr10_18cm_50MeV/CV080317_Zr10_18cm_50MeV_peak_data.csv'








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





df_AAZR01 = pd.read_csv(path_Zr+file_AAZR01)
df_ABZR09 = pd.read_csv(path_Zr+file_ABZR09)
df_ABZR08 = pd.read_csv(path_Zr+file_ABZR08)
df_ACZR07 = pd.read_csv(path_Zr+file_ACZR07)
df_AEZR06 = pd.read_csv(path_Zr+file_AEZR06)
df_AFZR10 = pd.read_csv(path_Zr+file_AFZR10)
df_AGZR09 = pd.read_csv(path_Zr+file_AGZR09)
df_AHZR08 = pd.read_csv(path_Zr+file_AHZR08)
df_AIZR07 = pd.read_csv(path_Zr+file_AIZR07)
df_AJZR06 = pd.read_csv(path_Zr+file_AJZR06)
df_AKZR06 = pd.read_csv(path_Zr+file_AKZR06)
df_ALZR06 = pd.read_csv(path_Zr+file_ALZR06)
df_AMZR06 = pd.read_csv(path_Zr+file_AMZR06)
df_ANZR07 = pd.read_csv(path_Zr+file_ANZR07)
df_AOZR07 = pd.read_csv(path_Zr+file_AOZR07)
df_APZR08 = pd.read_csv(path_Zr+file_APZR08)
df_AQZR09 = pd.read_csv(path_Zr+file_AQZR09)
df_ARZR10 = pd.read_csv(path_Zr+file_ARZR10)
df_ASZR07 = pd.read_csv(path_Zr+file_ASZR07)
df_AUZR08 = pd.read_csv(path_Zr+file_AUZR08)
df_AVZR09 = pd.read_csv(path_Zr+file_AVZR09)
df_AWZR10 = pd.read_csv(path_Zr+file_AWZR10)
df_BJZR06 = pd.read_csv(path_Zr+file_BJZR06)
df_BKZR07 = pd.read_csv(path_Zr+file_BKZR07)
df_BLZR08 = pd.read_csv(path_Zr+file_BLZR08)
df_BMZR09 = pd.read_csv(path_Zr+file_BMZR09)
df_BNZR10 = pd.read_csv(path_Zr+file_BNZR10)
df_BQZR06 = pd.read_csv(path_Zr+file_BQZR06)
df_CAZR06 = pd.read_csv(path_Zr+file_CAZR06)
df_CBZR07 = pd.read_csv(path_Zr+file_CBZR07)
df_CDZR08 = pd.read_csv(path_Zr+file_CDZR08)
df_CNZR09 = pd.read_csv(path_Zr+file_CNZR09)
df_CVZR10 = pd.read_csv(path_Zr+file_CVZR10)
df_CVZR10 = pd.read_csv(path_Zr+file_CVZR10)










df_concat_Zr = pd.concat((df_AXXR05, df_BAZR01, df_BIZR04, df_BPZR01, df_BRZR02, df_BSZR05, df_BTZR03, df_BXZR01, df_CHZR03, df_CRZR04, df_BFZR01, df_BGZR02, df_BHZR03, df_BOZR04, df_CYZR02, df_CZZR05, df_AAZR01, df_ABZR09, df_ABZR08, df_ACZR07, df_AEZR06, df_AFZR10, df_AGZR09, df_AHZR08, df_AIZR07, df_AJZR06, df_AKZR06, df_ALZR06, df_AMZR06, df_ANZR07, df_AOZR07, df_APZR08, df_AQZR09, df_ARZR10, df_ASZR07, df_AUZR08, df_AVZR09, df_AWZR10, df_BJZR06, df_BKZR07, df_BLZR08, df_BMZR09, df_BNZR10, df_BQZR06, df_CAZR06, df_CBZR07, df_CDZR08, df_CNZR09, df_CVZR10, df_CVZR10), axis = 0)
df_concat_Zr= df_concat_Zr[(df_concat_Zr['isotope'] != '90NB') | (df_concat_Zr['energy'] != 329.058)] # Exclude rows where isotope is '90NB' and energy is 329.058
df_concat_Zr.to_csv(path_Zr+file_concat_Zr)





def fit_prod_rate(isotope_list, isotope_chain_parent, foil, path, file, stack, plot=False):
    t_irr_h = 0.33

    R = {}
    for isotope in isotope_list:
        if isotope == '86ZR':
            R0 = 0
        else:
            R0 = 1E4
        R[isotope] =[[R0,t_irr_h]] 

    if stack == '30MeV':
        EoB = '02/13/2017 14:47:00'

    if stack == '50MeV':
        EoB = '02/12/2017 19:21:00'

    dc = ci.DecayChain(isotope_chain_parent, R=R, units='h')
    dc.get_counts(foil, EoB, path+file)
    isotopes, R, cov_R = dc.fit_R()
    print(f'A0 for {isotope_chain_parent}: ', dc.activity(isotope_chain_parent, 0))

    if plot==True:
        dc.plot()

    return isotopes, R, cov_R




def calc_prod_rates_in_foil(isotope_list_list, isotope_chain_parent_list, foil, stack, write_to_file=False, show_plot=False):
    print(' ')
    print(foil)
    csv_file_path = f'./Calculated_R/{foil}_R_by_curie.csv'
    if write_to_file==True:
            with open(csv_file_path, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerow(['Isotope', 'R', 'R_unc', 'Rel_R_unc_in_percent'])

    for isotope_list, parent in zip(isotope_list_list, isotope_chain_parent_list):
        isotopes, R, cov_R = fit_prod_rate(isotope_list, parent, foil, path_Zr, file_concat_Zr, stack, plot=show_plot)

        print(isotopes)
        print(R)
        # print(cov_R)
        if write_to_file==True:
            with open(csv_file_path, 'a', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                for i, isotope in enumerate(isotopes):
                    # print(i, isotope)
                    csv_writer.writerow([f'{isotope}', f'{R[i]}', f'{np.sqrt(cov_R[i,i])}', f'{np.sqrt(cov_R[i,i])/R[i]*100:.2e}'])




       

# isotope_list_list_Zr01 = [['87Y'], ['90NB'], ['96NB']]
# isotope_chain_parent_list_Zr01 = ['87Y', '90NB', '96NB']

# isotope_list_list_Zr02 = [['87Y'], ['90NB'], ['96NB']]
# isotope_chain_parent_list_Zr02 = ['87Y', '90NB', '96NB']

# isotope_list_list_Zr03 = [['90NB'], ['96NB']]
# isotope_chain_parent_list_Zr03 = ['90NB', '96NB']

# isotope_list_list_Zr04 = [['90NB'], ['96NB']]
# isotope_chain_parent_list_Zr04 = ['90NB', '96NB']

# isotope_list_list_Zr05 = [['90NB'], ['96NB']]
# isotope_chain_parent_list_Zr05 = ['90NB', '96NB']


# isotope_list_list_Zr06 = [['90NB'], ['96NB'],['86Y'], ['87Y']]
# isotope_chain_parent_list_Zr06 = ['90NB', '96NB', '86Y', '87Y']

# isotope_list_list_Zr07 = [['90NB'], ['96NB'],['86Y'], ['87Y']]
# isotope_chain_parent_list_Zr07 = ['90NB', '96NB', '86Y', '87Y']

# isotope_list_list_Zr08 = [['90NB'], ['96NB'], ['86Y'], ['87Y']]
# isotope_chain_parent_list_Zr08 = ['90NB', '96NB', '86Y', '87Y']

# isotope_list_list_Zr09 = [['90NB'], ['96NB'], ['86Y'], ['87Y']]
# isotope_chain_parent_list_Zr09 = ['90NB', '96NB', '86Y', '87Y']

# isotope_list_list_Zr10 = [['90NB'], ['96NB'], ['87Y']]
# isotope_chain_parent_list_Zr10 = ['90NB', '96NB', '87Y']



isotope_list_list_Zr06 = [['86Y']]
isotope_chain_parent_list_Zr06 = ['86Y']

isotope_list_list_Zr07 = [['86Y']]
isotope_chain_parent_list_Zr07 = ['86Y']

isotope_list_list_Zr08 = [['86Y']]
isotope_chain_parent_list_Zr08 = ['86Y']

isotope_list_list_Zr09 = [['86Y']]
isotope_chain_parent_list_Zr09 = ['86Y']





# isotope_list_list_Zr01 = [['87Ym', '87Y']]
# isotope_chain_parent_list_Zr01 = ['87Ym']

# isotope_list_list_Zr02 = [['87Ym', '87Y']]
# isotope_chain_parent_list_Zr02 = ['87Ym']

# isotope_list_list_Zr06 = [['87Ym', '87Y']]
# isotope_chain_parent_list_Zr06 = ['87Ym']

# isotope_list_list_Zr07 = [['87Ym', '87Y']]
# isotope_chain_parent_list_Zr07 = ['87Ym']

# isotope_list_list_Zr08 = [['87Ym', '87Y']]
# isotope_chain_parent_list_Zr08 = ['87Ym']

# isotope_list_list_Zr09 = [['87Ym', '87Y']]
# isotope_chain_parent_list_Zr09 = ['87Ym']

# isotope_list_list_Zr10 = [['87Ym', '87Y']]
# isotope_chain_parent_list_Zr10 = ['87Ym']





# isotope_list_list_Zr01 = [['88Y']]
# isotope_chain_parent_list_Zr01 = ['88Y']

# isotope_list_list_Zr02 = [['88Y']]
# isotope_chain_parent_list_Zr02 = ['88Y']

# isotope_list_list_Zr03 = [['88Y']]
# isotope_chain_parent_list_Zr03 = ['88Y']

# isotope_list_list_Zr04 = [['88Y']]
# isotope_chain_parent_list_Zr04 = ['88Y']

# isotope_list_list_Zr05 = [['88Y']]
# isotope_chain_parent_list_Zr05 = ['88Y']

# isotope_list_list_Zr06 = [['88NB', '88ZR', '88Y']]
# isotope_chain_parent_list_Zr06 = ['88NB']

# isotope_list_list_Zr07 = [['88NB', '88ZR', '88Y']]
# isotope_chain_parent_list_Zr07 = ['88NB']

# isotope_list_list_Zr08 = [['88NB', '88ZR', '88Y']]
# isotope_chain_parent_list_Zr08 = ['88NB']

# isotope_list_list_Zr09 = [['88NB', '88ZR', '88Y']]
# isotope_chain_parent_list_Zr09 = ['88NB']

# isotope_list_list_Zr10 = [['88NB', '88ZR', '88Y']]
# isotope_chain_parent_list_Zr10 = ['88NB']





# isotope_list_list_Zr01 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr01 = ['89NB']

# isotope_list_list_Zr02 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr02 = ['89NB']

# isotope_list_list_Zr03 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr03 = ['89NB']

# isotope_list_list_Zr04 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr04 = ['89NB']

# isotope_list_list_Zr05 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr05 = ['89NB']

# isotope_list_list_Zr06 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr06 = ['89NB']

# isotope_list_list_Zr07 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr07 = ['89NB']

# isotope_list_list_Zr08 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr08 = ['89NB']

# isotope_list_list_Zr09 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr09 = ['89NB']

# isotope_list_list_Zr10 = [['89NB', '89ZR']]
# isotope_chain_parent_list_Zr10 = ['89NB']



# isotope_list_list_Zr01 = [['90Ym']]
# isotope_chain_parent_list_Zr01 = ['90Ym']

# isotope_list_list_Zr03 = [['90Ym']]
# isotope_chain_parent_list_Zr03 = ['90Ym']

# isotope_list_list_Zr06 = [['90Ym']]
# isotope_chain_parent_list_Zr06 = ['90Ym']

# isotope_list_list_Zr07 = [['90Ym']]
# isotope_chain_parent_list_Zr07 = ['90Ym']

# isotope_list_list_Zr08 = [['90Ym']]
# isotope_chain_parent_list_Zr08 = ['90Ym']


# isotope_list_list_Zr10 = [['90Ym']]
# isotope_chain_parent_list_Zr10 = ['90Ym']





# isotope_list_list_Zr03 = [['91Y']]
# isotope_chain_parent_list_Zr03 = ['91Y']

# isotope_list_list_Zr04 = [['91Y']]
# isotope_chain_parent_list_Zr04 = ['91Y']

# isotope_list_list_Zr06 = [['91Y']]
# isotope_chain_parent_list_Zr06 = ['91Y']

# isotope_list_list_Zr07 = [['91Y']]
# isotope_chain_parent_list_Zr07 = ['91Y']

# isotope_list_list_Zr08 = [['91Y']]
# isotope_chain_parent_list_Zr08 = ['91Y']

# isotope_list_list_Zr09 = [['91Y']]
# isotope_chain_parent_list_Zr09 = ['91Y']

# isotope_list_list_Zr10 = [['91Y']]
# isotope_chain_parent_list_Zr10 = ['91Y']



# isotope_list_list_Zr01 = [['92Y']]
# isotope_chain_parent_list_Zr01 = ['92Y']

# isotope_list_list_Zr03 = [['92Y']]
# isotope_chain_parent_list_Zr03 = ['92Y']

# isotope_list_list_Zr06 = [['92Y']]
# isotope_chain_parent_list_Zr06 = ['92Y']

# isotope_list_list_Zr07 = [['92Y']]
# isotope_chain_parent_list_Zr07 = ['92Y']

# isotope_list_list_Zr08 = [['92Y']]
# isotope_chain_parent_list_Zr08 = ['92Y']



# isotope_list_list_Zr01 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr01 = ['95ZR']

# isotope_list_list_Zr02 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr02 = ['95ZR']

# isotope_list_list_Zr03 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr03 = ['95ZR']

# isotope_list_list_Zr04 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr04 = ['95ZR']

# isotope_list_list_Zr05 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr05 = ['95ZR']

# isotope_list_list_Zr06 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr06 = ['95ZR']

# isotope_list_list_Zr07 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr07 = ['95ZR']

# isotope_list_list_Zr08 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr08 = ['95ZR']

# isotope_list_list_Zr09 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr09 = ['95ZR']

# isotope_list_list_Zr10 = [['95ZR', '95NBm', '95NB']]
# isotope_chain_parent_list_Zr10 = ['95ZR']



# calc_prod_rates_in_foil(isotope_list_list_Zr01, isotope_chain_parent_list_Zr01, 'Zr01', '30MeV', write_to_file=True, show_plot=True)
# calc_prod_rates_in_foil(isotope_list_list_Zr02, isotope_chain_parent_list_Zr02, 'Zr02', '30MeV', write_to_file=True, show_plot=True)
# calc_prod_rates_in_foil(isotope_list_list_Zr03, isotope_chain_parent_list_Zr03, 'Zr03', '30MeV', write_to_file=True, show_plot=True)
# calc_prod_rates_in_foil(isotope_list_list_Zr04, isotope_chain_parent_list_Zr04, 'Zr04', '30MeV', write_to_file=True, show_plot=True)
# calc_prod_rates_in_foil(isotope_list_list_Zr05, isotope_chain_parent_list_Zr05, 'Zr05', '30MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Zr06, isotope_chain_parent_list_Zr06, 'Zr06', '50MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Zr07, isotope_chain_parent_list_Zr07, 'Zr07', '50MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Zr08, isotope_chain_parent_list_Zr08, 'Zr08', '50MeV', write_to_file=True, show_plot=True)
calc_prod_rates_in_foil(isotope_list_list_Zr09, isotope_chain_parent_list_Zr09, 'Zr09', '50MeV', write_to_file=True, show_plot=True)
# calc_prod_rates_in_foil(isotope_list_list_Zr10, isotope_chain_parent_list_Zr10, 'Zr10', '50MeV', write_to_file=True, show_plot=True)




